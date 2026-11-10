//! A round whose rail runs past the far end of the side edge it ends on, onto the next
//! face (SolidCraft's Knoten_laminiert Fillet4: a corner of a 0.5 mm plate rounded r 3
//! where the plate's top runs into a concave round 0.009 short of the rail).
//!
//! An L block: a plate 30 x 20 x 1 with a wall 6 high from x = 10, the concave edge
//! between the plate's top and the wall rounded r 2 (the plate's top runs into the round
//! tangentially at x = 8). The plate's vertical corner edge at the origin rounded R 8.01:
//! its rail on the front face (y = 0) is the line x = 8.01, past the end of the front
//! face's top edge with the plate's top (x = 8) by 0.01, where the r 2 round has risen
//! 2 - sqrt(4 - 0.01^2). The rail ends on the next edge, the front face's with the r 2
//! round, as Fusion ends it.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::boolean::api_brep::boolean_union;
use solidcore::curve::curve::CurveKind;
use solidcore::kernel::UpstreamFixes;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::brep::p1;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::topology::{edge, vertex, Brep, EdgeId};
use solidcore::Kernel;

fn edge_between(b: &Brep, a: Point3, c: Point3) -> EdgeId {
    let at = |p: Point3, q: Point3| p.distance(q) < 1e-6;
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

#[test]
fn a_rail_past_its_side_edge_ends_on_the_next_one() {
    let mut k = Kernel::default();
    // SolidCraft's set: every fix but the two it holds back
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k.fixes.fillet_round_ends_on_outline = false;
    let plate = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 30.0, 20.0, 1.0).unwrap();
    let wall = create_box(&mut k, Point3::new(10.0, 0.0, 0.0), 20.0, 20.0, 6.0).unwrap();
    let mut b = boolean_union(&mut k, plate, wall).unwrap();
    let valley = edge_between(
        &b,
        Point3::new(10.0, 0.0, 1.0),
        Point3::new(10.0, 20.0, 1.0),
    );
    fillet_edges(&mut k, &mut b, &[valley], 2.0, 1, 1, 1.0).unwrap();
    let corner = edge_between(&b, Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    let r = 8.01;
    fillet_edges(&mut k, &mut b, &[corner], r, 1, 1, 1.0).unwrap();

    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert!(check.is_valid(), "{:?}", check.lines());
    // the rail's end: on the front face's edge with the r 2 round
    let end = Point3::new(r, 0.0, 3.0 - (4.0 - (r - 8.0) * (r - 8.0)).sqrt());
    let near = p1::vertices(&b, None)
        .unwrap()
        .into_iter()
        .map(|v| vertex::point(&b, v).unwrap().distance(end))
        .fold(f64::INFINITY, f64::min);
    assert!(near < 1e-6, "{near}");
}

/// Where the round's end crosses the plate's top edge with the r 2 round, 6.2e-6 from that
/// edge's end (both corners, the stitch gluing the round's vertex to the edge's end from either
/// side): the vertex is the crossing, on the round, the plate's top and the r 2 round, as Fusion
/// has it, not the edge's old end (the vertices lie within each other's tolerance).
#[test]
fn the_rounds_end_crosses_the_seam_where_it_lies() {
    let mut k = Kernel::default();
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k.fixes.fillet_round_ends_on_outline = false;
    let plate = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 30.0, 20.0, 1.0).unwrap();
    let wall = create_box(&mut k, Point3::new(10.0, 0.0, 0.0), 20.0, 20.0, 6.0).unwrap();
    let mut b = boolean_union(&mut k, plate, wall).unwrap();
    let valley = edge_between(
        &b,
        Point3::new(10.0, 0.0, 1.0),
        Point3::new(10.0, 20.0, 1.0),
    );
    fillet_edges(&mut k, &mut b, &[valley], 2.0, 1, 1, 1.0).unwrap();
    let r = 8.01;
    let corners = [
        edge_between(&b, Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0)),
        edge_between(&b, Point3::new(0.0, 20.0, 0.0), Point3::new(0.0, 20.0, 1.0)),
    ];
    fillet_edges(&mut k, &mut b, &corners, r, 1, 1, 1.0).unwrap();
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert!(check.is_valid(), "{:?}", check.lines());
    let points: Vec<Point3> = p1::vertices(&b, None)
        .unwrap()
        .into_iter()
        .map(|v| vertex::point(&b, v).unwrap())
        .collect();
    let nearest = |p: Point3| {
        points
            .iter()
            .map(|q| q.distance(p))
            .fold(f64::INFINITY, f64::min)
    };
    let dy = r - (r * r - (r - 8.0) * (r - 8.0)).sqrt();
    for (crossing, old_end) in [
        (Point3::new(8.0, dy, 1.0), Point3::new(8.0, 0.0, 1.0)),
        (
            Point3::new(8.0, 20.0 - dy, 1.0),
            Point3::new(8.0, 20.0, 1.0),
        ),
    ] {
        assert!(
            nearest(crossing) < 1e-9,
            "{crossing:?}: {}",
            nearest(crossing)
        );
        assert!(nearest(old_end) > 1e-9, "{old_end:?}: {}", nearest(old_end));
    }
}

/// Every edge ends on its vertices, those at the crossings too: put there, a vertex takes the
/// edges that ended at the seam's old end with it (Knoten_laminiert's straight seam edges ended 1.4e-5 past their vertex,
/// and a sketch projecting one drew it so, where Fusion's line ends on the vertex).
#[test]
fn the_edges_at_the_crossing_end_on_it() {
    let mut k = Kernel::default();
    k.fixes = UpstreamFixes::all();
    k.fixes.fillet_refuse_rounds_breaking_through = false;
    k.fixes.fillet_round_ends_on_outline = false;
    let plate = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 30.0, 20.0, 1.0).unwrap();
    let wall = create_box(&mut k, Point3::new(10.0, 0.0, 0.0), 20.0, 20.0, 6.0).unwrap();
    let mut b = boolean_union(&mut k, plate, wall).unwrap();
    let valley = edge_between(
        &b,
        Point3::new(10.0, 0.0, 1.0),
        Point3::new(10.0, 20.0, 1.0),
    );
    fillet_edges(&mut k, &mut b, &[valley], 2.0, 1, 1, 1.0).unwrap();
    let r = 8.01;
    let corners = [
        edge_between(&b, Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0)),
        edge_between(&b, Point3::new(0.0, 20.0, 0.0), Point3::new(0.0, 20.0, 1.0)),
    ];
    fillet_edges(&mut k, &mut b, &corners, r, 1, 1, 1.0).unwrap();
    // Every edge of the body ends on its vertices: a line on the vertex, a fitted curve (the
    // round's) within 1e-8.
    let mut off = Vec::new();
    for e in p1::edges(&b, None).unwrap() {
        let curve = edge::curve(&b, e).unwrap().unwrap();
        let line = matches!(curve.kind(), CurveKind::Line(_));
        let ivl = edge::interval(&b, e).unwrap();
        let ends = [
            (edge::start_vertex(&b, e).unwrap(), ivl.min()),
            (edge::end_vertex(&b, e).unwrap(), ivl.max()),
        ];
        for (v, t) in ends {
            let at = vertex::point(&b, v.unwrap()).unwrap();
            let d = curve.evaluate_point(&k, t).unwrap().distance(at);
            if d >= if line { 1e-9 } else { 1e-8 } {
                off.push(format!("{e:?} (line {line}) {d:e} from {at:?}"));
            }
        }
    }
    assert!(off.is_empty(), "{}", off.join("\n"));
}
