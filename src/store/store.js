import { create } from 'zustand';
import { devtools } from 'zustand/middleware'; // 1. devtools를 import 합니다.

const parameter = create(
  // 2. 기존 스토어 로직을 devtools로 감싸줍니다.
  devtools(
    (set) => ({
      nodes: [], // 빈 배열로 초기화
      members: [],
      sections: [],
      materials: [],
      forces: [],
      displacements: [],

      // 상태 변경 함수도 필요하면 추가
      setNodes: (nodes) => set({ nodes }),
      setMembers: (members) => set({ members }),
      setSections: (sections) => set({ sections }),
      setMaterials: (materials) => set({ materials }),
      setForce: (forces) => set({ forces }),
      setDisplacement: (displacements) => set({ displacements }),
    }),
    { name: 'ParameterStore' } // 3. DevTools에 표시될 스토어 이름을 지정합니다.
  )
);

// 전역에 할당 (개발 편의를 위해)
if (typeof window !== 'undefined') {
  window.__STORE__ = parameter;
}

export default parameter;