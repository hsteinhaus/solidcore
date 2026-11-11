//! Not upstream (`fixes.fillet_blend_vertices_on_their_faces`): a round's vertices lie on the
//! faces they bound, and their edges end on them.
//!
//! SolidCraft's Knoten_laminiert after its Shell1 (`data/knoten_laminiert_shell1.smb`): two
//! tubes crossing on a 0.5 mm plate, rounded r 2 at their roots, shelled 0.5 inward. Its
//! Fillet4 rounds the plate outline's eight vertical corner edges at the arm ends r 3 (45
//! degree arcs, 0.5 high). Each round's top end meets the tube root's round where that runs
//! tangentially into the plate top. Made at the ends of the rounds' fitted rails, those
//! vertices lie 2.09e-5 off the root round, and SMLib's check refuses the body (`SmVertex`
//! test #1, `Bad Vertex Gap`, against 2e-5). With the fix every vertex lies on its faces and
//! the body passes the check at level 2.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::geom::attribute::AttributeStore;
use solidcore::kernel::{Kernel, UpstreamFixes};
use solidcore::math::Point3;
use solidcore::query::api_queries::edge_closest_point;
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::{loopuse, EdgeId};
use solidcore::Brep;

const SHELL1: &[u8] = include_bytes!("data/knoten_laminiert_shell1.smb");

/// Fillet4's picks: the middles of the eight corner edges.
const CORNERS: [[f64; 3]; 8] = [
    [-45.608387386532335, -25.102290732122416, 0.25],
    [-25.10229073212245, -45.60838738653231, 0.25],
    [25.102290732122437, -45.60838738653233, 0.25],
    [45.60838738653233, -25.102290732122455, 0.25],
    [45.608387386532335, 25.102290732122416, 0.25],
    [25.102290732122444, 45.60838738653231, 0.25],
    [-45.60838738653233, 25.10229073212245, 0.25],
    [-25.102290732122434, 45.608387386532314, 0.25],
];

/// Every opt-in fix on.
fn kernel() -> Kernel {
    Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    }
}

/// A SolidCraft file, read under the notification naming as SolidCraft reads it.
fn read(k: &mut Kernel, bytes: &[u8]) -> Brep {
    let mut attributes = AttributeStore::default();
    let mut records = NotifyRecords::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &loopuse::notify_ref,
        records: &mut records,
        next_key: 1 << 50,
    };
    let mut b = read_brep_from_file(k, Some(bytes.to_vec()), false, false, &mut env).unwrap();
    b.refresh_cached_topology();
    b
}

/// The edge nearest `p`, within 0.1.
fn pick(k: &Kernel, b: &Brep, p: [f64; 3]) -> EdgeId {
    let q = Point3::new(p[0], p[1], p[2]);
    let mut best = (0.1f64, None);
    for &e in b.edge_list() {
        if let Ok(c) = edge_closest_point(k, b, e, q) {
            if c.closest.distance <= best.0 {
                best = (c.closest.distance, Some(e));
            }
        }
    }
    best.1.unwrap_or_else(|| panic!("no edge near {p:?}"))
}

#[test]
fn knoten_fillet4_vertices_lie_on_their_faces() {
    let mut k = kernel();
    let mut b = read(&mut k, SHELL1);
    let edges: Vec<EdgeId> = CORNERS.iter().map(|&p| pick(&k, &b, p)).collect();
    fillet_edges(&mut k, &mut b, &edges, 3.0, 1, 1, 0.0).unwrap();
    assert_eq!(
        (
            b.face_list().len(),
            b.edge_list().len(),
            b.vertex_list().len()
        ),
        (64, 152, 90),
        "Fusion's Fillet4"
    );
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
}
