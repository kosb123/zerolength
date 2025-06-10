import { create } from 'zustand';

const parameter = create(set => ({
  nodes: [],         // 빈 배열로 초기화
  members: [],
  sections: [],
  materials: [],
  forces: [],
  displacements: [],
  variables: [],

  // 상태 변경 함수도 필요하면 추가
  setNodes: (nodes) => set({ nodes }),
  setMembers: (members) => set({ members }),
  setSections: (sections) => set({ sections }),
  setMaterials: (materials) => set({ materials }),
  setForce: (forces) => set({ forces }),
  setDisplacement: (displacements) => set({ displacements }),
  setVariables: (variables) => set({ variables }),
  addVariable: (variable) => set(state => ({ variables: [...state.variables, variable] })),
}));

// 전역에 할당 (개발 편의를 위해)
if (typeof window !== 'undefined') {
  window.__STORE__ = parameter;
}

export default parameter;
