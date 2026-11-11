//! Not upstream (opt-in `Kernel::fixes.fillet_blend_vertices_on_their_faces`): the vertices
//! of a fillet's blend faces put where the surfaces of their faces meet.
//!
//! A vertex the fillet makes is the end of a fitted rail or section curve, so it lies off a
//! face it bounds by up to that curve's fit ([`put_vertices_on_their_faces`]).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletSession::put_blend_vertices_on_their_faces`] (crate) | none (own design) |

use super::finish::FinishAccess;
use super::FilletSession;
use crate::error::{EntityRef, Result};
use crate::geom::attribute::ids;
use crate::kernel::Kernel;
use crate::topology::brep::p1;
use crate::topology::brep::vertices_on_faces::put_vertices_on_their_faces;
use crate::topology::{vertex, FaceId, VertexId};

impl FilletSession<'_> {
    /// [`put_vertices_on_their_faces`] for the vertices of the blend faces (marked
    /// [`ids::FILLET_BLEND_FACE`]).
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
            put_vertices_on_their_faces(k, b, env, &vertices)
        })
    }
}
