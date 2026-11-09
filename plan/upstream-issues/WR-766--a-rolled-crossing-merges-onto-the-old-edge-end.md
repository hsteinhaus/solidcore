#### WR-766. A round's end crossing a seam is merged onto the seam's old end, off the crossing by its tolerance
- **Severity:** Wrong result (geometry within tolerance, but not the exact point: Knoten_laminiert Fillet4 seam vertices 1.4e-5 off Fusion's; later sketch projections onto them differ)
- **Area:** Fillet finish (merge of the fillet into the body)
- **Location:** Rust: `blend/fillet_executive/session/finish.rs` (`place_rolled_crossings`, after `finish_merge`).
- **What happens:** the vertex WR-765 makes where a round's end crosses the end face's edge with the next face lies within tolerance of that edge's old end. The merge keeps whichever of the two points its order meets first, so on Knoten 8 of 16 such vertices stay at the old end (1.4e-5 off).
- **Fix:** after the merge, a body vertex within its tolerance of such a fillet-built crossing is put on the crossing when that point lies nearer all the vertex's faces than its own.
- **Confirmed:** `tests/fillet_rail_end_past_side_edge.rs::the_rounds_end_crosses_the_seam_where_it_lies`; Knoten_laminiert Fillet4's 16 seam vertices at Fusion's to 8e-14.
