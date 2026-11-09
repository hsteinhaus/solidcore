#### WR-763. Deleting a wire whose ends come free makes its faceuse shell a vertex shell
- **Severity:** Crash (a dangling shell: `NotFound` on `topology::Brep::shell`)
- **Area:** Topology (`SmBrep::DeleteEdge`, its wire branch)
- **Location:** `SmBrep::DeleteEdge`: with both ends of a wire edge shell vertices afterwards (or its one end, for a closed wire), the edge's shell is reset to the vertex shell of the first end. Rust: `crates/solidcore/src/topology/brep/p4.rs`, `delete_wire_edge`.
- **What happens:** the reset takes the shell for a wire shell holding the edge alone. A faceuse shell can own a wire whose ends touch only other shells' wires: the lamina edges of faces `MakeManifold` deletes, once the wires at their ends are deleted before them (SolidCraft's Knoten_laminiert Shell1, after the merge). Made a vertex shell, the faceuse shell lists none of its face uses, which keep naming it; `MakeManifold` then deletes the vertex with its shell, and the next query of a face use's shell fails.
- **Fix:** only a wire shell is reset; a freed end of a faceuse shell's wire gets a vertex shell of its own, as the branch with one free end already does.
- **Confirmed:** `topology::brep::p4::tests::a_free_wire_of_a_faceuse_shell_leaves_the_shell_its_faces` (a cube's outside shell owning a free wire: the shell keeps its six face uses, each end a vertex shell of its own).
