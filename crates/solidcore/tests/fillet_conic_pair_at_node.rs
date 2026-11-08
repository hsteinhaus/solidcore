//! A round whose two quadrics meet in a pair of plane conics, classified at a node of the
//! pair (where both conics pass), follows the conic of its own edge (WR-757). SolidCraft's
//! Knoten_laminiert (`tests/data/solidcraft_knoten_laminiert_fillet2.smb`, the body as
//! SolidCraft wrote it): four half-tubes r 11.5 on a 0.5 mm plate, along the diagonals,
//! crossing at the origin; each concave valley between two arms runs from the crossing's
//! top `(0, 0, 11.5)` down to the arms' root round (Fillet1) and on over it to the plate.
//! Rounded r 2 (Fusion's Fillet2), each valley adds a round and its ball's end, 5.1874 mm3.
//! The two arms' offsets meet in the ellipses of the planes `x = 0` and `y = 0`, and the
//! round is classified at its first point, the crossing's node, on both: the plane was the
//! first one tried, the other valley's for one of the four, and that round failed (status
//! 1001, a rollover read off rails built round the wrong ellipse).
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
use solidcore::kernel::UpstreamFixes;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::{p1, p10::TopologyEnv, p2};
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::{edge, vertex, Brep, EdgeId};
use solidcore::{EntityRef, Kernel};

const BODY: &[u8] = include_bytes!("data/solidcraft_knoten_laminiert_fillet2.smb");

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

/// The edge from `a` to `b` (either way).
fn edge_between(b: &Brep, a: Point3, c: Point3) -> EdgeId {
    let at = |p: Point3, q: Point3| p.distance(q) < 1e-3;
    let found: Vec<EdgeId> = p1::edges(b, None)
        .unwrap()
        .into_iter()
        .filter(|&e| {
            let s = vertex::point(b, edge::start_vertex(b, e).unwrap().unwrap()).unwrap();
            let t = vertex::point(b, edge::end_vertex(b, e).unwrap().unwrap()).unwrap();
            (at(s, a) && at(t, c)) || (at(s, c) && at(t, a))
        })
        .collect();
    assert_eq!(found.len(), 1, "{a:?} {c:?}");
    found[0]
}

/// The valley toward `(x, y)` (a unit axis direction in the plane): the arms' crossing
/// ellipse and its run on over the root round to the plate.
fn valley(b: &Brep, x: f64, y: f64) -> [EdgeId; 2] {
    let top = Point3::new(0.0, 0.0, 11.5);
    let root = Point3::new(15.982 * x, 15.982 * y, 2.129_63);
    let plate = Point3::new(18.762 * x, 18.762 * y, 0.5);
    [edge_between(b, top, root), edge_between(b, root, plate)]
}

/// The volume and faces rounding the valleys toward `dirs` r 2 adds.
fn round(dirs: &[(f64, f64)]) -> (f64, usize) {
    let mut k = Kernel::default();
    // SolidCraft's set: every fix but the two it holds back
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k.fixes.fillet_round_ends_on_outline = false;
    let mut b = body(&mut k);
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let faces = p1::faces(&b, None).unwrap().len();
    let before = brep_volume(&k, &mut b, 1e-10, &attrs).unwrap();
    let edges: Vec<EdgeId> = dirs.iter().flat_map(|&(x, y)| valley(&b, x, y)).collect();
    fillet_edges(&mut k, &mut b, &edges, 2.0, 1, 1, 1.0)
        .unwrap_or_else(|e| panic!("{dirs:?}: {e:?}"));
    assert!(p2::is_manifold_solid(&b).unwrap(), "{dirs:?}");
    let added = brep_volume(&k, &mut b, 1e-10, &attrs).unwrap() - before;
    (added, p1::faces(&b, None).unwrap().len() - faces)
}

const VALLEYS: [(f64, f64); 4] = [(0.0, 1.0), (1.0, 0.0), (0.0, -1.0), (-1.0, 0.0)];

#[test]
fn each_valley_s_round_follows_its_own_conic() {
    for dir in VALLEYS {
        let (added, faces) = round(&[dir]);
        assert_eq!(faces, 2, "{dir:?}");
        assert!((added - 5.1874).abs() < 1e-4, "{dir:?}: {added}");
    }
}

/// Fusion's Fillet2 whole, the four valleys in one call (WR-758): at the crossing's top every
/// round has lost its width, the arms' faces tangent there, and the rounds meet in the point
/// with no corner patch. The adjacent valleys' corner there was refused (an Nx2 corner, no
/// topology). 25 + 8 faces, 39383.603 + 4 x 5.1874 = 39404.352 mm3, Fusion's.
#[test]
fn the_four_valleys_round_as_fusion_s_fillet2() {
    let (added, faces) = round(&VALLEYS);
    assert_eq!(faces, 8);
    assert!((added - 4.0 * 5.1874).abs() < 4e-4, "{added}");
}
