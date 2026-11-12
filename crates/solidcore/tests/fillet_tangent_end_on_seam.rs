//! Not upstream (`fixes.fillet_trace_points_in_their_edges_period` in the fork's
//! `sm_TangentEndRib`): a round's end rib on a periodic end face runs within one period.
//!
//! SolidCraft's Pixracer_Basis before its Fillet2 (`data/pixracer_basis_before_fillet2.smb`):
//! R1.5 on the four straight inner edges of the walls' top, each running on tangentially into
//! a corner cylinder (r 7, full turn, its seam on the wall's tangent line). Each round ends on
//! that cylinder in a rib from the tangent end (on the seam, v = 0) to the other rail (v =
//! 42.43, one period from -1.55): the rib's first guess crossed the whole cylinder and no
//! station solved (1001). With the fix the four rounds are built and the body passes the check
//! at level 2.
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

const INPUT: &[u8] = include_bytes!("data/pixracer_basis_before_fillet2.smb");

/// Fillet2's picks: the middles of the four inner top edges.
const EDGES: [[f64; 3]; 4] = [
    [0.0, 35.0, 36.8],
    [-35.0, 0.0, 36.8],
    [0.0, -35.0, 36.8],
    [35.0, 0.0, 36.8],
];

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

/// The edge nearest `p`, within 0.01.
fn pick(k: &Kernel, b: &Brep, p: [f64; 3]) -> EdgeId {
    let q = Point3::new(p[0], p[1], p[2]);
    let mut best = (0.01, None);
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
fn pixracer_fillet2_ends_on_the_corner_cylinders() {
    let mut k = Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    };
    let mut b = read(&mut k, INPUT);
    let edges: Vec<EdgeId> = EDGES.iter().map(|&p| pick(&k, &b, p)).collect();
    fillet_edges(&mut k, &mut b, &edges, 1.5, 1, 1, 0.0).unwrap();
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
}
