import { useRef, useEffect } from 'react';
import { useFrame } from '@react-three/fiber';
import { OrbitControls, Text, Billboard } from '@react-three/drei';
import * as THREE from 'three';
import parameter from './store/store';

function MyElement3D() {
  const refMesh = useRef();

  // Zustand store에서 상태 읽기
  const nodes = parameter(state => state.nodes || []);
  const members = parameter(state => state.members || []);

  // 노드ID별 좌표를 빠르게 찾기 위한 Map 생성
  const nodeMap = new Map();
  nodes.forEach(({ id, x, y, z }) => {
    nodeMap.set(id, new THREE.Vector3(x, y, z));
  });

  const controlsRef = useRef();

  // Alt 키 누름 여부에 따라 OrbitControls 모드 변경
  useEffect(() => {
    // MOUSE ACTION CONSTANTS: ROTATE=0, DOLLY=1, PAN=2
    const MOUSE_ROTATE = 0;
    const MOUSE_PAN = 2;

    const updateControls = (isAltPressed) => {
      if (controlsRef.current) {
        // Alt가 눌려있으면 회전(0), 아니면 이동(2)
        const newAction = isAltPressed ? MOUSE_ROTATE : MOUSE_PAN;

        // 상태가 바뀔 때만 업데이트
        if (controlsRef.current.mouseButtons.MIDDLE !== newAction) {
          controlsRef.current.mouseButtons.MIDDLE = newAction;
          controlsRef.current.update();
        }
      }
    };

    const handleKeyChange = (e) => {
      updateControls(e.altKey);
    };

    // 포커스를 잃었을 때(Alt 탭 등) 상태 초기화
    const handleBlur = () => {
      updateControls(false);
    };

    window.addEventListener('keydown', handleKeyChange);
    window.addEventListener('keyup', handleKeyChange);
    window.addEventListener('blur', handleBlur);

    // 초기 상태 강제 설정 (Pan)
    updateControls(false);

    return () => {
      window.removeEventListener('keydown', handleKeyChange);
      window.removeEventListener('keyup', handleKeyChange);
      window.removeEventListener('blur', handleBlur);
    };
  }, []);

  useFrame(() => {
    if (refMesh.current) {
      // 애니메이션 처리 가능
    }
  });

  return (
    <>
      <directionalLight position={[5, 5, 5]} intensity={1} />
      <axesHelper args={[1000]} />
      <OrbitControls
        ref={controlsRef}
        enableDamping={true}
        dampingFactor={1}
        mouseButtons={{
          LEFT: null,
          MIDDLE: 2, // 2 = PAN
          RIGHT: null
        }}
        listenToKeyEvents={window}
      />

      {/* 노드 표시 (구 + 번호) */}
      {nodes.map(({ id, x, y, z }, idx) => (
        <group key={idx} position={[x, y, z]}>
          <mesh>
            <sphereGeometry args={[0.05, 16, 16]} />
            <meshBasicMaterial color="black" />
          </mesh>
          <Billboard>
            <Text
              position={[0.15, 0.15, 0]}
              fontSize={0.15}
              color="black"
              anchorX="left"
              anchorY="bottom"
            >
              {id}
            </Text>
          </Billboard>
        </group>
      ))}

      {/* 멤버별 연결선 그리기 */}
      {members.map((elem, idx) => {
        const start = nodeMap.get(elem.n1);
        const end = nodeMap.get(elem.n2);

        if (!start || !end) return null;

        const points = [start, end];
        const geometry = new THREE.BufferGeometry().setFromPoints(points);

        return (
          <line key={idx} geometry={geometry}>
            <lineBasicMaterial attach="material" color="orange" linewidth={2} />
          </line>
        );
      })}
    </>
  );
}

export default MyElement3D;
