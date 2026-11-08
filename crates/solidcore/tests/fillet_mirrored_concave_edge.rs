//! A concave round ending at a mixed-convexity Nx1 corner whose end face meets more than
//! the round and one other face is built the same on a body and on its mirror image
//! (WR-755). SolidCraft's Motortraeger (`tests/data/solidcraft_motortraeger_mirrored_fillet.smb`,
//! the body as SolidCraft wrote it): a wall `y = 36` and two struts, one modelled at
//! `x = 19` and its mirror image at `x = 1` (about `x = 10`), meeting it in concave edges
//! from `z = 0` to `z = 11`. At the bottom end both side edges are cliff edges, and the end
//! face the corner extends was the first side's, which is the bottom face on one side and
//! the strut's sloped face on the other: the bottom face does not hold the other side
//! edge, the corner fell back to a Coons patch, and with the patch keeping the fillet's
//! edges (`fillet_patch_keeps_fillet_edges`) the stitch left the body open.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::geom::vector3d::Point3;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::{p1, p10::TopologyEnv, p2};
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::{edge, vertex, Brep, EdgeId};
use solidcore::{EntityRef, Kernel};

const BODY: &[u8] = include_bytes!("data/solidcraft_motortraeger_mirrored_fillet.smb");

fn body(k: &mut Kernel) -> Brep {
    let key = |_: EntityRef| ObjectRef {
        key: ObjectKey(0),
        type_code: 0,
    };
    let mut attributes = AttributeStore::default();
    let mut records = NotifyRecords::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &key,
        records: &mut records,
        next_key: 1 << 40,
    };
    let mut b = read_brep_from_file(k, Some(BODY.to_vec()), false, false, &mut env).unwrap();
    b.use_arena_naming().unwrap();
    b.refresh_cached_topology();
    b
}

/// The concave line edge from `(x, 36, 0)` to `(x, 36, 11)`.
fn edge_at(b: &Brep, x: f64) -> EdgeId {
    let at = |p: Point3, z: f64| p.distance(Point3::new(x, 36.0, z)) < 1e-6;
    let found: Vec<EdgeId> = p1::edges(b, None)
        .unwrap()
        .into_iter()
        .filter(|&e| {
            let s = vertex::point(b, edge::start_vertex(b, e).unwrap().unwrap()).unwrap();
            let t = vertex::point(b, edge::end_vertex(b, e).unwrap().unwrap()).unwrap();
            (at(s, 0.0) && at(t, 11.0)) || (at(s, 11.0) && at(t, 0.0))
        })
        .collect();
    assert_eq!(found.len(), 1, "x {x}");
    found[0]
}

/// The volume and faces a round of 3 on the edge at `x` adds.
fn round(x: f64) -> (f64, usize) {
    let mut k = Kernel::default();
    k.fixes.fillet_patch_keeps_fillet_edges = true;
    let mut b = body(&mut k);
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let faces = p1::faces(&b, None).unwrap().len();
    let before = brep_volume(&k, &mut b, 1e-10, &attrs).unwrap();
    let e = edge_at(&b, x);
    fillet_edges(&mut k, &mut b, &[e], 3.0, 1, 1, 1.0).unwrap_or_else(|e| panic!("x {x}: {e:?}"));
    assert!(p2::is_manifold_solid(&b).unwrap(), "x {x}");
    let added = brep_volume(&k, &mut b, 1e-10, &attrs).unwrap() - before;
    (added, p1::faces(&b, None).unwrap().len() - faces)
}

#[test]
fn a_mirrored_concave_round_is_built_as_its_original() {
    let (original, original_faces) = round(19.0);
    let (mirrored, mirrored_faces) = round(1.0);
    assert_eq!((original_faces, mirrored_faces), (1, 1));
    assert!((original - 20.827).abs() < 1e-3, "original {original}");
    assert!(
        (mirrored - original).abs() < 1e-9,
        "mirrored {mirrored} original {original}"
    );
}
