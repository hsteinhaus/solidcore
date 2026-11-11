//! Not upstream (opt-in `Kernel::fixes.fillet_blend_vertices_on_their_faces`): the vertices
//! of a fillet's blend faces put where the surfaces of their faces meet.
//!
//! A vertex the fillet makes is the end of a fitted rail or section curve, so it lies off a
//! face it bounds by up to that curve's fit: off by more than the B-rep's intersection
//! tolerance, `SmVertex` test #1 (`Bad Vertex Gap`) refuses the body. Where three faces meet
//! at the vertex, their three surfaces meet in one point near it, found from the vertex's
//! drops onto them; the vertex moves there when that point is nearer all three surfaces
//! than the vertex is, and its edges' ends follow it.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletSession::put_blend_vertices_on_their_faces`] (crate) | none (own design) |

use super::finish::FinishAccess;
use super::FilletSession;
use crate::error::{EntityRef, Result};
use crate::geom::attribute::ids;
use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::SolverOperation;
use crate::geom::tol::{eff_zero, ANG_TOL_RAD};
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::Vec3;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::nurbs::knots::Side;
use crate::surface::surface::p4::NoObject;
use crate::surface::surface::p7::base_evaluate_normal;
use crate::surface::surface::Surface;
use crate::topology::brep::p1;
use crate::topology::edge::NotifyEnv;
use crate::topology::{face, vertex, vertexuse, Brep, FaceId, VertexId};

impl FilletSession<'_> {
    /// Puts each vertex of the blend faces (marked [`ids::FILLET_BLEND_FACE`]) where three
    /// faces meet onto their surfaces' common point ([`put_on_its_faces`]), and the ends of
    /// every blend face vertex's edges onto the vertex
    /// ([`p1::update_vertex_edge_intervals`]).
    ///
    /// # Errors
    ///
    /// Those of the topology reads, the drops, the solve, the move and the notifications.
    pub(super) fn put_blend_vertices_on_their_faces(&mut self, k: &mut Kernel) -> Result<()> {
        let mut blends: Vec<FaceId> = Vec::new();
        for f in p1::faces(self.objects.target, None)? {
            let key = self.target_identity(EntityRef::Face(f)).key;
            if self
                .objects
                .attributes
                .find_attribute(key, ids::FILLET_BLEND_FACE)
                .is_some()
            {
                blends.push(f);
            }
        }
        let mut vertices: Vec<VertexId> = Vec::new();
        for v in p1::vertices(self.objects.target, None)? {
            if vertex::faces(self.objects.target, v)?
                .iter()
                .any(|f| blends.contains(f))
            {
                vertices.push(v);
            }
        }
        FinishAccess(&mut self.objects).with_target_env(k, |k, b, env| {
            for v in vertices {
                put_on_its_faces(k, b, env, v)?;
                // its edges end on it, whether it moved or not (a fitted rail's or section's
                // end is off the vertex it was made at by up to the fit)
                p1::update_vertex_edge_intervals(k, b, env, v)?;
                for vu in vertex::vertexuses(b, v)? {
                    vertexuse::init_vertex_face_gap(b, vu)?;
                }
            }
            Ok(())
        })
    }
}

/// Moves `v` to the point nearest the surfaces of its three faces ([`nearest_to_all`] from
/// its drops onto them, their meeting point where they cross), when that is nearer all
/// three than `v` is and within `v`'s tolerance of it. A vertex of other than three faces, or
/// one the solve does not improve, stays where it is.
fn put_on_its_faces(k: &Kernel, b: &mut Brep, env: &mut NotifyEnv<'_>, v: VertexId) -> Result<()> {
    let mut faces = vertex::faces(b, v)?;
    faces.sort_unstable();
    faces.dedup();
    let [f1, f2, f3] = faces[..] else {
        return Ok(());
    };
    let at = vertex::point(b, v)?;
    let mut surfaces: Vec<Surface> = Vec::new();
    let mut domains = Vec::new();
    let mut guesses: Vec<Point2> = Vec::new();
    let mut gap = 0.0_f64;
    let uses = vertex::vertexuses(b, v)?;
    for f in [f1, f2, f3] {
        // the vertex's drop onto the face as its use there makes it (`SmVertex` test #1's)
        let mut drop = None;
        for &vu in &uses {
            if let Some(fu) = vertexuse::faceuse(b, vu)? {
                if b.faceuse(fu)?.face == Some(f) {
                    drop = Some(vertexuse::compute_uv_point(k, b, vu, false)?);
                    break;
                }
            }
        }
        let (Some(uv), Some(surface)) = (drop, face::p2::surface(b, f)?) else {
            return Ok(());
        };
        let surface = surface.clone_for_evaluation()?;
        // expanded as `FilletRemoval::CalcNewCornerVertices` expands its srf/srf/srf
        // domains: the meeting point of a vertex at a face's boundary may lie past it
        let mut domain = surface.natural_uv_domain_with_nurb(k)?;
        domain.expand_relative(2.0);
        gap = gap.max(surface.evaluate_point(k, uv)?.distance(at));
        surfaces.push(surface);
        domains.push(domain);
        guesses.push(uv);
    }
    let Some((to, spread)) = nearest_to_all(k, &surfaces, &domains, &mut guesses, at)? else {
        return Ok(());
    };
    if !(spread < gap) || to.distance(at) > vertex::tolerance(b, v)? {
        return Ok(());
    }
    let key_of = env.key_of;
    vertex::set_point(
        b,
        v,
        to,
        env.records,
        env.attributes,
        env.context.as_deref_mut(),
        key_of(EntityRef::Vertex(v)),
        &|vu| key_of(EntityRef::Vertexuse(vu)),
    )?;
    Ok(())
}

/// The most Gauss-Newton steps of [`nearest_to_all`] (it converges quadratically where the
/// faces cross; a step count past this is a solve going nowhere).
const MAX_STEPS: usize = 32;

/// The point nearest all of `surfaces` near `from` (least squares of its distances to them),
/// with its largest distance to any: Gauss-Newton steps, each the least move onto the
/// tangent planes at the point's drops onto the surfaces (from `uvs`, updated), two planes
/// within the angle tolerance of parallel (faces tangent there) taken as one. Where the
/// surfaces cross transversally that is the point they meet in; where two are tangent, the
/// point on the third nearest both. `None` when a drop finds nothing.
fn nearest_to_all(
    k: &Kernel,
    surfaces: &[Surface],
    domains: &[Extent2d],
    uvs: &mut [Point2],
    from: Point3,
) -> Result<Option<(Point3, f64)>> {
    let mut x = from;
    let mut spread = f64::INFINITY;
    for _ in 0..MAX_STEPS {
        // the drops and their normals
        let mut planes: Vec<(Vec3, f64)> = Vec::new();
        spread = 0.0;
        for ((surface, domain), uv) in surfaces.iter().zip(domains).zip(uvs.iter_mut()) {
            let r = surface.local_point_solve::<NoObject>(
                k,
                domain,
                SolverOperation::Minimize,
                x,
                *uv,
            )?;
            if !r.found {
                return Ok(None);
            }
            let start = &r.solution.start;
            *uv = Point2::new(start.param(0)?, start.param(1)?);
            let p = surface.evaluate_point(k, *uv)?;
            spread = spread.max(p.distance(x));
            let Ok(n) = base_evaluate_normal(k, surface, *uv, Side::Left, Side::Left)?.unitize()
            else {
                return Ok(None);
            };
            planes.push((n, n.dot(p - x)));
        }
        // the least move d onto the planes: Gram-Schmidt on their normals, a normal within
        // the angle tolerance of the earlier ones' span dropped
        let mut basis: Vec<(Vec3, f64)> = Vec::new();
        for (n, r) in planes {
            let mut q = n;
            let mut c = r;
            for &(b, cb) in &basis {
                let a = q.dot(b);
                q = q - b * a;
                c -= a * cb;
            }
            let len = q.length();
            if len <= ANG_TOL_RAD.sin() {
                continue;
            }
            basis.push((q * (1.0 / len), c / len));
        }
        let d = basis
            .iter()
            .fold(Vec3::new(0.0, 0.0, 0.0), |d, &(q, c)| d + q * c);
        x = x + d;
        if d.length() <= eff_zero() * (1.0 + x.max_dimension()) {
            break;
        }
    }
    Ok(Some((x, spread)))
}
