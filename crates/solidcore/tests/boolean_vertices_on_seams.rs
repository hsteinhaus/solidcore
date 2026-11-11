//! Not upstream (`fixes.boolean_vertices_on_seams`): a Boolean's vertices within their
//! tolerance of a periodic face's seam lie on it.
//!
//! SolidCraft's Pixracer_Basis at its Combine1 (`data/pixracer_combine1_{a,b}.smb`): a plate
//! with four corner cones (r 7, seams at x = ±35 / y = ±28 and turned) joined with a shelled
//! body whose flat walls x = ±35 / y = ±35 touch the cones tangentially along their seams.
//! Where a wall's edge meets a seam the union leaves the vertex 1e-7 off it, and SMLib's check
//! refuses the four cone faces (`SmFace` test #36, `Face Near missing Seam`). With the fix
//! each such vertex lies on its seam and the body passes the check at level 2.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::boolean::api_brep::boolean;
use solidcore::geom::a_object::BooleanOperation;
use solidcore::geom::attribute::AttributeStore;
use solidcore::kernel::{Kernel, UpstreamFixes};
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::loopuse;
use solidcore::topology::topology::NotifyRecords;
use solidcore::Brep;

const A: &[u8] = include_bytes!("data/pixracer_combine1_a.smb");
const B: &[u8] = include_bytes!("data/pixracer_combine1_b.smb");

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

#[test]
fn pixracer_combine1_vertices_lie_on_the_cone_seams() {
    let mut k = kernel();
    let a = read(&mut k, A);
    let b = read(&mut k, B);
    let mut r = boolean(&mut k, a, b, BooleanOperation::Union)
        .unwrap()
        .result;
    let lines = check_brep(&mut k, &mut r, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
}
