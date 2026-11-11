//! Not upstream: vertices put where the surfaces of their faces meet, and their edges' ends
//! onto them.
//!
//! An operation that makes a vertex at the end of a fitted curve (a fillet's rail, an offset's
//! trimmed rim) leaves it off a face it bounds by up to that curve's fit, and its edges ending
//! off it by as much: off by more than the B-rep's intersection tolerance `SmVertex` test #1
//! (`Bad Vertex Gap`) refuses the body; within it, a later operation can take the vertex for
//! one it only nearly meets (`SmFace` test #36, a seam's near miss). Where three faces meet at
//! the vertex, their three surfaces meet in one point near it, found from the vertex's drops
//! onto them; the vertex moves there when that point is nearer all three surfaces than the
//! vertex is, and its edges' ends follow it.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`put_vertices_on_their_faces`] | none (own design) |
//! | [`put_vertices_on_seams`] | none (own design; upstream's `Fix_BadGaps` for `m_eBadNearMissSeam` was never built) |

use crate::error::{EntityRef, Result};
use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::SolverOperation;
use crate::geom::tol::{eff_zero, ANG_TOL_RAD};
use crate::geom::vector2d::Point2;
use crate::geom::vector3d::Vec3;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::nurbs::b_spline_surface::SurfParamType;
use crate::nurbs::knots::Side;
use crate::surface::surface::p4::NoObject;
use crate::surface::surface::p7::base_evaluate_normal;
use crate::surface::surface::Surface;
use crate::topology::brep::p1;
use crate::topology::edge::NotifyEnv;
use crate::topology::{face, vertex, vertexuse, Brep, VertexId};

/// Puts each of `vertices` where three faces meet onto their surfaces' common point
/// ([`put_on_its_faces`]), and the ends of each one's edges onto it
/// ([`p1::update_vertex_edge_intervals`]), whether it moved or not.
///
/// # Errors
///
/// Those of the topology reads, the drops, the move, the edge updates and the notifications.
pub fn put_vertices_on_their_faces(
    k: &Kernel,
    b: &mut Brep,
    env: &mut NotifyEnv<'_>,
    vertices: &[VertexId],
) -> Result<()> {
    for &v in vertices {
        put_on_its_faces(k, b, env, v)?;
        // its edges end on it, whether it moved or not (a fitted curve's end is off the vertex
        // it was made at by up to the fit)
        p1::update_vertex_edge_intervals(k, b, env, v)?;
        for vu in vertex::vertexuses(b, v)? {
            vertexuse::init_vertex_face_gap(b, vu)?;
        }
    }
    Ok(())
}

/// Puts each of `vertices` that lies within its tolerance of the seam of a periodic face it
/// bounds exactly onto that seam ([`put_on_its_seam`]), and the ends of its edges onto it
/// ([`p1::update_vertex_edge_intervals`]). Two operands' boundaries that meet where a face of
/// one touches the other's periodic face tangentially along its seam (a plane tangent to a
/// cone along the generator its seam lies on) meet there only within the tolerance, and the
/// vertex is left off the seam: `SmFace` test #36 (`Face Near missing Seam`) refuses it.
///
/// # Errors
///
/// Those of the topology reads, the drops, the move, the edge updates and the notifications.
pub fn put_vertices_on_seams(
    k: &Kernel,
    b: &mut Brep,
    env: &mut NotifyEnv<'_>,
    vertices: &[VertexId],
) -> Result<()> {
    for &v in vertices {
        if put_on_its_seam(k, b, env, v)? {
            p1::update_vertex_edge_intervals(k, b, env, v)?;
            for vu in vertex::vertexuses(b, v)? {
                vertexuse::init_vertex_face_gap(b, vu)?;
            }
        }
    }
    Ok(())
}

/// Moves `v` onto the nearest seam of a periodic face it bounds (its drop there with the seam
/// parameter put in), when that is within `v`'s tolerance of it and keeps it within that of
/// every face it bounds. Whether it moved.
fn put_on_its_seam(k: &Kernel, b: &mut Brep, env: &mut NotifyEnv<'_>, v: VertexId) -> Result<bool> {
    let at = vertex::point(b, v)?;
    let tol = vertex::tolerance(b, v)?;
    // its drops onto its faces
    let mut drops: Vec<(Surface, Extent2d, Point2)> = Vec::new();
    for vu in vertex::vertexuses(b, v)? {
        let Some(fu) = vertexuse::faceuse(b, vu)? else {
            continue;
        };
        let Some(f) = b.faceuse(fu)?.face else {
            continue;
        };
        let Some(surface) = face::p2::surface(b, f)? else {
            continue;
        };
        let surface = surface.clone_for_evaluation()?;
        let domain = surface.natural_uv_domain_with_nurb(k)?;
        let uv = vertexuse::compute_uv_point(k, b, vu, false)?;
        drops.push((surface, domain, uv));
    }
    // the nearest seam point within its tolerance
    let mut to = None;
    let mut nearest = tol;
    for (surface, domain, uv) in &drops {
        for dir in [SurfParamType::U, SurfParamType::V] {
            if !surface.is_closed(k, domain, dir, None, false)?.0 {
                continue;
            }
            let (lo, hi, t) = if dir == SurfParamType::U {
                (domain.min().x, domain.max().x, uv.x)
            } else {
                (domain.min().y, domain.max().y, uv.y)
            };
            let seam = if t - lo <= hi - t { lo } else { hi };
            let on = if dir == SurfParamType::U {
                Point2::new(seam, uv.y)
            } else {
                Point2::new(uv.x, seam)
            };
            let p = surface.evaluate_point(k, on)?;
            let d = p.distance(at);
            if d > 0.0 && d <= nearest {
                nearest = d;
                to = Some(p);
            }
        }
    }
    let Some(to) = to else {
        return Ok(false);
    };
    // still within its tolerance of every face
    for (surface, domain, uv) in &drops {
        let mut domain = *domain;
        domain.expand_relative(2.0);
        let r = surface.local_point_solve::<NoObject>(
            k,
            &domain,
            SolverOperation::Minimize,
            to,
            *uv,
        )?;
        let start = &r.solution.start;
        if !r.found
            || surface
                .evaluate_point(k, Point2::new(start.param(0)?, start.param(1)?))?
                .distance(to)
                > tol
        {
            return Ok(false);
        }
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
    Ok(true)
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
