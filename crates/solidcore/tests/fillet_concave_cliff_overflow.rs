//! Not upstream (`fixes.fillet_concave_cliff_overflow`): a concave round whose rails leave
//! the face it rolls on across cliffs adds its region cut at the cliffs.
//!
//! SolidCraft's Pixracer_Basis before its Fillet3 (`data/pixracer_basis_before_fillet3.smb`):
//! R1.5 where four corner towers (r 8.2, tangent to the outer walls) stand on a plate at z 38.
//! Each round runs 270 degrees round its tower; at both ends its rail on the plate leaves it
//! past the walls x, y = -/+36.2, and on the way it crosses a hole in the plate. Fusion cuts the
//! round flush with the walls and the hole's wall. Upstream refuses it (`ConcaveOverflow`).
//! With the fix the four regions are united: 7 faces more, 69.24 mm3 more (Fusion: 69.25), and
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
use solidcore::mass::integrator::properties::brep_volume_default;
use solidcore::math::Point3;
use solidcore::query::api_queries::edge_closest_point;
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::{loopuse, EdgeId};
use solidcore::Brep;

const INPUT: &[u8] = include_bytes!("data/pixracer_basis_before_fillet3.smb");

/// Fillet3's picks: the middles of the four tower roots.
const EDGES: [[f64; 3]; 4] = [
    [22.202, 22.202, 38.0],
    [22.202, -22.202, 38.0],
    [-22.202, -22.202, 38.0],
    [-22.202, 22.202, 38.0],
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
fn pixracer_fillet3_rounds_the_towers_flush_with_the_walls() {
    let mut k = Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    };
    let mut b = read(&mut k, INPUT);
    let edges: Vec<EdgeId> = EDGES.iter().map(|&p| pick(&k, &b, p)).collect();
    let faces = b.face_list().len();
    let before = brep_volume_default(&k, &mut b, 1e-6).unwrap();
    fillet_edges(&mut k, &mut b, &edges, 1.5, 1, 1, 0.0).unwrap();
    assert_eq!(b.face_list().len(), faces + 7);
    let added = brep_volume_default(&k, &mut b, 1e-6).unwrap() - before;
    assert!((added - 69.248).abs() < 0.01, "added {added}");
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
}
