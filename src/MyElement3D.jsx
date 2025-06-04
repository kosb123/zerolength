import { useRef } from 'react';
import { useFrame } from '@react-three/fiber';
import { OrbitControls, Text } from '@react-three/drei';
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

  useFrame(() => {
    if (refMesh.current) {
      // 애니메이션 처리 가능
    }
  });

  return (
    <>
      <directionalLight position={[5, 5, 5]} intensity={1} />
      <axesHelper args={[1000]} />
      <OrbitControls enableDamping={true} dampingFactor={1} />

      {/* 노드 표시 (구 + 번호) */}
      {nodes.map(({ id, x, y, z }, idx) => (
        <group key={idx} position={[x, y, z]}>
          <mesh>
            <sphereGeometry args={[0.1, 16, 16]} />
            <meshBasicMaterial color="skyblue" />
          </mesh>
          <Text
            position={[0.15, 0, 0]}
            fontSize={0.15}
            color="black"
            anchorX="left"
            anchorY="middle"
          >
            {id}
          </Text>
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
