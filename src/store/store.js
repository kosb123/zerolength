import { create } from 'zustand';

const useStore = create((set) => ({
  nodes: [],
  members: [],
  materials: [], // Combined Material & Section Properties

  // Boundary Conditions & Loads (Solver Input)
  constraints: [], // Array of constrained DOFs (1-based), e.g., [1, 2, 3, ...]
  forces: [],      // Flat force vector (size = 6 * nodeCount)

  // Analysis Results
  analysisResults: {
    displacements: null, // Flat vector
    reactions: null,     // Flat vector
    elementForces: [],   // Array of { elem_id, f_local }
  },

  // Actions
  setNodes: (nodes) => set({ nodes }),
  setMembers: (members) => set({ members }),
  setMaterials: (materials) => set({ materials }),
  setConstraints: (constraints) => set({ constraints }),
  setForces: (forces) => set({ forces }),

  setAnalysisResults: (results) => set({ analysisResults: results }),

  reset: () => set({
    nodes: [], members: [], materials: [],
    constraints: [], forces: [],
    analysisResults: { displacements: null, reactions: null, elementForces: [] }
  }),
}));

// Expose store to window for debugging
if (typeof window !== 'undefined') {
  window.__STORE__ = useStore;
}

export default useStore;
