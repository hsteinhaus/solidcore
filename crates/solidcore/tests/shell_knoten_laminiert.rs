//! SolidCraft's Knoten_laminiert (`tests/data/solidcraft_knoten_laminiert_shell1.smb`, the
//! body as SolidCraft wrote it): four half-tubes along the diagonals of a 0.5 mm plate,
//! crossing at the origin, their four concave valleys rounded r 2 (Fillet2), each round
//! ending at the crossing's top in a pole. Fusion's Shell1 opens the bottom and the four
//! arm ends and shells the rest 0.5 inward.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::UpstreamFixes;
use solidcore::offset::offset_executive::p2::{shell_brep_simple, OffsetOptions, ShellOptions};
use solidcore::topology::api_files::read_brep_from_file;
use solidcore::topology::brep::{p1, p10::TopologyEnv};
use solidcore::topology::topology::NotifyRecords;
use solidcore::topology::Brep;
use solidcore::{EntityRef, Kernel, SolidError};

const BODY: &[u8] = include_bytes!("data/solidcraft_knoten_laminiert_shell1.smb");

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

/// Shell1 as SolidCraft runs it: its fixes (every one but the two it holds back), the
/// bottom and the four arm ends open (SolidCraft's picks: faces 17, 0, 15, 7, 5), 0.5 inward.
fn shell1() -> Result<Brep, SolidError> {
    let mut k = Kernel::default();
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k.fixes.fillet_round_ends_on_outline = false;
    let b = body(&mut k);
    let faces = p1::faces(&b, None).unwrap();
    let open: Vec<_> = [17usize, 0, 15, 7, 5].iter().map(|&i| faces[i]).collect();
    let options = ShellOptions {
        offset: OffsetOptions {
            distance: -0.5,
            extend_convex_edges: true,
            merge_results: true,
        },
        create_offset_solid: true,
        merge_sheet_side_walls: true,
    };
    shell_brep_simple(&mut k, &b, options, &open)
}

/// The rounds' offset faces keep their UV trim curves when put on their extended surfaces
/// (WR-760): the extension is the offset surface over its domain, so the curves the face
/// has are its curves there. Dropped anew, the rails' curves went astray at the pole their
/// rails start from (any v of the collapsed side is that point: v -81.8), their trim
/// domain ran past the extension's, and the shell failed in its extended offset
/// (`ExtendResultFaces`, status 1001). The shell now gets past that step (its merge still
/// fails at the node, a separate defect).
#[test]
fn the_rounds_offsets_keep_their_trims_on_their_extensions() {
    let made = shell1();
    if let Err(SolidError::Upstream { operation, .. }) = &made {
        assert_ne!(
            *operation, "offset::offset_executive::p2::do_extended_offset",
            "{made:?}"
        );
    }
}

/// Shell1 builds as Fusion's: 56 faces, 4603.714 mm³, a manifold solid (with WR-759 to
/// WR-763: the pole made exact, the trims kept on the extensions, the offset of the pole
/// evaluated, the slits removed, the free wires' shells kept).
#[test]
fn shell1_is_fusions() {
    use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
    let mut out = shell1().unwrap();
    assert_eq!(p1::faces(&out, None).unwrap().len(), 56);
    assert!(solidcore::topology::brep::p2::is_manifold_solid(&out).unwrap());
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let volume = brep_volume(&Kernel::default(), &mut out, 1e-8, &attrs).unwrap();
    assert!((volume - 4603.714).abs() < 1e-3, "{volume}");
}
