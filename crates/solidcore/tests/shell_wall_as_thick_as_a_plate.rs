//! A shell whose wall is as thick as a plate it opens onto (WR-762): a 20 x 20 x 0.5 plate
//! with a 10 x 10 x 10 tower on it, its bottom open, shelled 0.5 inward. The plate's top,
//! offset by the wall, lands on the open bottom's plane, and the offsets of the plate's
//! four sides keep no height but at their ends. The tower is hollowed through the plate;
//! the plate stays as it is: 20 x 20 x 0.5 + 10 x 10 x 10 - 9 x 9 x 10 = 390, 16 faces (the
//! plate's 4 sides, its top and bottom rings, the tower's 4 outer and 4 inner walls and its
//! top's two faces). As SolidCraft's Knoten_laminiert's Shell1 does at its 0.5 plate.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::boolean::api_brep::{boolean_difference, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::offset::offset_executive::p2::{shell_brep_simple, OffsetOptions, ShellOptions};
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::brep::{p1, p2};
use solidcore::Kernel;

/// The plate (`plate` thick) and tower shelled 0.5 inward, their bottom open, with
/// SolidCraft's fixes (every one but the two it holds back): the result's faces and
/// volume.
fn shell(plate: f64) -> (usize, f64) {
    shell_cut(plate, false)
}

/// [`shell`], with `chamfer`: the plate's corner at (20, 20) cut off at 45 degrees, 3 along
/// each side.
fn shell_cut(plate: f64, chamfer: bool) -> (usize, f64) {
    let mut k = Kernel::default();
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k.fixes.fillet_round_ends_on_outline = false;
    let mut base = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 20.0, 20.0, plate).unwrap();
    if chamfer {
        // a square of diagonal 6 about the corner, turned 45 degrees
        let side = 6.0 / 2f64.sqrt();
        let mut cut = create_box(
            &mut k,
            Point3::new(-side / 2.0, -side / 2.0, -1.0),
            side,
            side,
            plate + 2.0,
        )
        .unwrap();
        let (c, s) = (
            std::f64::consts::FRAC_PI_4.cos(),
            std::f64::consts::FRAC_PI_4.sin(),
        );
        let turn = [
            [c, -s, 0.0, 20.0],
            [s, c, 0.0, 20.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        solidcore::topology::brep::mirror::transform_matrix(&mut k, &mut cut, &turn, false)
            .unwrap();
        base = boolean_difference(&mut k, base, cut).unwrap();
    }
    let tower = create_box(&mut k, Point3::new(5.0, 5.0, plate), 10.0, 10.0, 10.0).unwrap();
    let b = boolean_union(&mut k, base, tower).unwrap();
    // the bottom, at z 0
    let bottom: Vec<_> = p1::faces(&b, None)
        .unwrap()
        .into_iter()
        .filter(|&f| {
            let s = b.face(f).unwrap().surface.as_ref().unwrap();
            let d = solidcore::topology::face::p2::uv_domain(&b, f).unwrap();
            s.evaluate_point(&k, d.mid()).unwrap().z.abs() < 1e-9
        })
        .collect();
    assert_eq!(bottom.len(), 1);
    let options = ShellOptions {
        offset: OffsetOptions {
            distance: -0.5,
            extend_convex_edges: true,
            merge_results: true,
        },
        create_offset_solid: true,
        merge_sheet_side_walls: true,
    };
    let mut made = shell_brep_simple(&mut k, &b, options, &bottom).unwrap();
    assert!(p2::is_manifold_solid(&made).unwrap());
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let v = brep_volume(&k, &mut made, 1e-8, &attrs).unwrap();
    (p1::faces(&made, None).unwrap().len(), v)
}

#[test]
fn a_wall_as_thick_as_the_plate_leaves_the_plate() {
    let (faces, v) = shell(0.5);
    assert_eq!(faces, 16);
    assert!((v - 390.0).abs() < 1e-6, "{v}");
}

/// A plate thicker than the wall keeps a cavity under its top, 19 x 19 x 0.5: no offset
/// face has a slit, and the shell is as it was: 20 x 20 x 1 + 10 x 10 x 10 - 19 x 19 x 0.5
/// - 9 x 9 x 10 = 409.5, 21 faces.
#[test]
fn a_plate_thicker_than_the_wall_is_shelled_as_before() {
    let (faces, v) = shell(1.0);
    assert_eq!(faces, 21);
    assert!((v - 409.5).abs() < 1e-6, "{v}");
}

/// The same with the plate's corner cut off at 45 degrees: there the cut's offset runs out
/// along the bottom past where the next side's offset comes back along it (a slit along
/// part of an edge). 20 x 20 - 4.5 = 395.5 of plate, 0.5 thick, + 1000 - 810 = 387.75,
/// 17 faces.
#[test]
fn a_wall_as_thick_as_a_cut_plate_leaves_the_plate() {
    let (faces, v) = shell_cut(0.5, true);
    assert_eq!(faces, 17);
    assert!((v - 387.75).abs() < 1e-6, "{v}");
}
