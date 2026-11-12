//! Pixracer_Basis before its Fillet1 (the body SolidCraft passed to the refused blend, binary SMB):
//! the R0.8 round between the floor (z = 1.2) and the rib (an R12.7 cylinder on the
//! diagonal), ending on the wall x = -35, is a four-sided fillet face `SmBrep::DeleteFilletFace`
//! removes: the floor and the rib then meet in their sharp edge, which reaches the wall where
//! the plane meets the cylinder there, (-35, -35 + 17.799, 1.2) (|x - y| = sqrt(12.7² - 1.7²) √2).
//! The input passes SMLib's check (level 2), and so does the result. Needs
//! `fixes.fillet_removal_continues_analytic_faces`: the rib and the floor are revolved surfaces
//! (cones) that only continue as themselves past their profiles' ends.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::blend::fillet_removal::delete_fillet_face;
use solidcore::geom::attribute::AttributeStore;
use solidcore::kernel::{Kernel, UpstreamFixes};
use solidcore::math::Point3;
use solidcore::query::api_queries::{brep_is_manifold_solid, find_edge, get_edges};
use solidcore::surface::plane::Plane;
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::face::p2 as face2;
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::{edge, loopuse, vertex};
use solidcore::Brep;

const INPUT: &[u8] = include_bytes!("data/pixracer_basis_before_fillet1.smb");

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
fn the_round_between_floor_and_rib_is_removed_to_their_sharp_edge() {
    let mut k = Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    };
    let mut b = read(&mut k, INPUT);
    let faces = b.face_list().len();
    assert!(check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines()
        .is_empty());
    // The wall's edge on the round, and the round beside it.
    let pick = find_edge(
        &k,
        &b,
        Point3::new(-35.0, -16.960_473_707_785, 1.384_159_712_864),
    )
    .unwrap();
    let arc = get_edges(&b).unwrap()[pick.index];
    let round = edge::faces(&b, arc, None)
        .unwrap()
        .into_iter()
        .find(|&f| {
            !face2::surface(&b, f)
                .unwrap()
                .unwrap()
                .is_kind_of(Plane::TYPE)
                .unwrap()
        })
        .unwrap();
    b.set_editing_enabled(true);
    delete_fillet_face(&mut k, &mut b, round, None).unwrap();
    assert_eq!(b.face_list().len(), faces - 1);
    assert!(brep_is_manifold_solid(&b).unwrap());
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
    let corner = Point3::new(
        -35.0,
        -35.0 + (12.7f64 * 12.7 - 1.7 * 1.7).sqrt() * 2f64.sqrt(),
        1.2,
    );
    let near = b.vertices().ids().any(|v| {
        let p = vertex::point(&b, v).unwrap();
        ((p.x - corner.x).powi(2) + (p.y - corner.y).powi(2) + (p.z - corner.z).powi(2)).sqrt()
            < 1e-6
    });
    assert!(
        near,
        "no vertex where floor and rib meet the wall, {corner:?}"
    );
}

/// All eight R0.8 rounds between the floors and the corner cones come off, one after another:
/// each corner cone grows twice, once for the round on either side, and stays a cone.
#[test]
fn every_round_between_floor_and_corner_comes_off_one_after_another() {
    use solidcore::blend::fillet_removal::FilletRemoval;
    use solidcore::blend::identify_fillet::{identify_fillet, IdentifyFilletOptions};
    let mut k = Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    };
    let mut b = read(&mut k, INPUT);
    let faces = b.face_list().len();
    let mut rounds = Vec::new();
    for f in b.face_list().to_vec() {
        let options = IdentifyFilletOptions {
            tol_3d: face2::tolerance(&b, f).unwrap(),
            tangency_tol_deg: 2.0,
            require_tangent: true,
            require_manifold: true,
            want_min_max_radii: false,
        };
        let id = identify_fillet(&k, &mut b, f, &options).unwrap();
        if id.is_fillet && id.solver_type == 0 && (id.radius - 0.8).abs() < 1e-6 {
            rounds.push(f);
        }
    }
    assert_eq!(rounds.len(), 8);
    b.set_editing_enabled(true);
    for f in rounds {
        let mut removal = FilletRemoval::new(f, face2::tolerance(&b, f).unwrap());
        removal.remove_face(&mut k, &mut b).unwrap();
        b.refresh_cached_topology();
        let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
            .unwrap()
            .lines();
        assert!(lines.is_empty(), "after {f:?}: {lines:?}");
    }
    assert_eq!(b.face_list().len(), faces - 8);
    assert!(brep_is_manifold_solid(&b).unwrap());
}

/// Fillet1's R1.5 edges once the R0.8 rounds are off: the wall bottoms, the walls' rib edges and
/// the corner rims (midpoints).
const FILLET1: [[f64; 3]; 16] = [
    [-35.0, -20.72753230511315, 7.209625992822557],
    [-35.0, 0.0, 1.2],
    [-35.0, 20.72753230823282, 7.209625995710183],
    [-32.94966009728894, 32.94983483774862, 12.199999999381042],
    [-20.727532305113165, 35.0, 7.209625992822565],
    [0.0, 35.0, 1.2],
    [20.72753230823279, 35.0, 7.209625995710187],
    [32.94967809521779, 32.949816840373835, 12.199999999594302],
    [35.0, 20.727532305113158, 7.2096259928225415],
    [35.0, 0.0, 1.2],
    [35.0, -20.727532308232817, 7.209625995710216],
    [32.94966009728892, -32.94983483774863, 12.199999999381035],
    [20.727532305113147, -35.0, 7.20962599282255],
    [0.0, -35.0, 1.2],
    [-20.727532308177445, -35.0, 7.2096259956589375],
    [-32.949834837748625, -32.94966009728891, 12.199999999381033],
];

/// Not upstream (`fillet_trace_points_in_their_edges_period`, the plane UVs of
/// `surface_plane_curve` about the plane's own domain, the blend vertices on their seams): with
/// the eight rounds off, Fillet1's R1.5 on all sixteen edges. A corner rim ends on its corner
/// cylinder's seam, and its trace started a period off (1001); a rib round's corner curve on the
/// floor, a plane whose domain starts at (84, 37), got its UVs about (0, 0), 3.85 off the floor
/// (the merge refused it). The body passes the check at level 2.
#[test]
fn fillet1_rounds_the_walls_once_the_rounds_are_off() {
    use solidcore::blend::fillet_executive::api::fillet_edges;
    use solidcore::blend::fillet_removal::FilletRemoval;
    use solidcore::blend::identify_fillet::{identify_fillet, IdentifyFilletOptions};
    use solidcore::query::api_queries::edge_closest_point;
    let mut k = Kernel {
        fixes: UpstreamFixes::all(),
        ..Kernel::default()
    };
    let mut b = read(&mut k, INPUT);
    let mut rounds = Vec::new();
    for f in b.face_list().to_vec() {
        let options = IdentifyFilletOptions {
            tol_3d: face2::tolerance(&b, f).unwrap(),
            tangency_tol_deg: 2.0,
            require_tangent: true,
            require_manifold: true,
            want_min_max_radii: false,
        };
        let id = identify_fillet(&k, &mut b, f, &options).unwrap();
        if id.is_fillet && id.solver_type == 0 && (id.radius - 0.8).abs() < 1e-6 {
            rounds.push(f);
        }
    }
    b.set_editing_enabled(true);
    for f in rounds {
        let mut removal = FilletRemoval::new(f, face2::tolerance(&b, f).unwrap());
        removal.remove_face(&mut k, &mut b).unwrap();
        b.refresh_cached_topology();
    }
    let edges: Vec<_> = FILLET1
        .iter()
        .map(|p| {
            let q = Point3::new(p[0], p[1], p[2]);
            let mut best = (0.01, None);
            for &e in b.edge_list() {
                if let Ok(c) = edge_closest_point(&k, &b, e, q) {
                    if c.closest.distance <= best.0 {
                        best = (c.closest.distance, Some(e));
                    }
                }
            }
            best.1.unwrap_or_else(|| panic!("no edge near {p:?}"))
        })
        .collect();
    fillet_edges(&mut k, &mut b, &edges, 1.5, 1, 1, 0.0).unwrap();
    assert_eq!(b.face_list().len(), 76);
    let lines = check_brep(&mut k, &mut b, CheckLevel::Level2)
        .unwrap()
        .lines();
    assert!(lines.is_empty(), "{lines:?}");
}
