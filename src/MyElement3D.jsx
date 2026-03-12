import { useRef, useEffect, useMemo } from 'react';
import { useFrame } from '@react-three/fiber';
import { OrbitControls, Text, Billboard } from '@react-three/drei';
import * as THREE from 'three';
import parameter from './store/store';

// 단일 Draw Call로 구체(노드)를 렌더링하여 성능을 대폭 높이는 InstancedMesh 컴포넌트
function NodeInstancedMesh({ nodes, color, size }) {
  const meshRef = useRef();
  const dummy = useMemo(() => new THREE.Object3D(), []);

  useFrame(({ camera }) => {
    if (!meshRef.current || nodes.length === 0) return;
    
    nodes.forEach((node, i) => {
      dummy.position.set(node.x, node.y, node.z);
      
      const distance = camera.position.distanceTo(dummy.position);
      const scale = distance * 0.01 * size;
      
      dummy.scale.set(scale, scale, scale);
      dummy.updateMatrix();
      meshRef.current.setMatrixAt(i, dummy.matrix);
    });
    meshRef.current.instanceMatrix.needsUpdate = true;
  });

  if (nodes.length === 0) return null;

  return (
    <instancedMesh ref={meshRef} args={[null, null, nodes.length]}>
      <sphereGeometry args={[1, 16, 16]} />
      <meshBasicMaterial color={color} />
    </instancedMesh>
  );
}

// 텍스트(ID) 표시는 노드마다 내용이 달라 개별 렌더링합니다.
function NodeLabel({ pos, id, color, size }) {
  const groupRef = useRef();

  useFrame(({ camera }) => {
    if (groupRef.current) {
      const distance = camera.position.distanceTo(groupRef.current.position);
      const scale = distance * 0.01 * size;
      groupRef.current.scale.set(scale, scale, scale);
    }
  });

  return (
    <group ref={groupRef} position={pos}>
      <Billboard>
        <Text position={[1.5, 1.5, 0]} fontSize={1.5} color={color} anchorX="left" anchorY="bottom">
          {id}
        </Text>
      </Billboard>
    </group>
  );
}

function MyElement3D() {
  const refMesh = useRef();

  // Zustand store에서 상태 읽기
  const nodes = parameter(state => state.nodes || []);
  const members = parameter(state => state.members || []);
  const nodeSettings = parameter(state => state.nodeSettings || { size: 1, color: '#000000' });
  
  const analysisResults = parameter(state => state.analysisResults);
  const deformationScale = parameter(state => state.deformationScale);

  // 노드ID별 좌표를 빠르게 찾기 위한 Map 생성 (useMemo로 최적화)
  const nodeMap = useMemo(() => {
    const map = new Map();
    nodes.forEach(({ id, x, y, z }) => {
      map.set(id, new THREE.Vector3(x, y, z));
    });
    return map;
  }, [nodes]);

  // 여러 개의 선을 한 번의 Draw Call로 그리기 위해 geometry 최적화 (Original Shape)
  const linesGeometry = useMemo(() => {
    if (members.length === 0 || nodeMap.size === 0) return null;
    const positions = [];
    members.forEach((elem) => {
      const start = nodeMap.get(elem.n1);
      const end = nodeMap.get(elem.n2);
      if (start && end) {
        positions.push(start.x, start.y, start.z);
        positions.push(end.x, end.y, end.z);
      }
    });

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    return geometry;
  }, [members, nodeMap]);

  // 변형된 노드 좌표 계산 (Deformed Shape)
  const { deformedNodes, deformedLinesGeometry } = useMemo(() => {
    if (!analysisResults || !analysisResults.displacements || nodes.length === 0) {
      return { deformedNodes: [], deformedLinesGeometry: null };
    }

    const { displacements } = analysisResults;
    const deformedMap = new Map();
    const dNodes = [];

    nodes.forEach(node => {
      // 0-based index for DOFs. u_x, u_y, u_z are the first 3 DOFs for a node.
      const idx = (node.id - 1) * 6;
      const ux = displacements[idx + 0] || 0;
      const uy = displacements[idx + 1] || 0;
      const uz = displacements[idx + 2] || 0;

      const dx = node.x + ux * deformationScale;
      const dy = node.y + uy * deformationScale;
      const dz = node.z + uz * deformationScale;

      dNodes.push({ id: node.id, x: dx, y: dy, z: dz });
      deformedMap.set(node.id, new THREE.Vector3(dx, dy, dz));
    });

    const positions = [];
    members.forEach((elem) => {
      const start = deformedMap.get(elem.n1);
      const end = deformedMap.get(elem.n2);
      if (start && end) {
        positions.push(start.x, start.y, start.z);
        positions.push(end.x, end.y, end.z);
      }
    });

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));

    return { deformedNodes: dNodes, deformedLinesGeometry: geometry };
  }, [nodes, members, analysisResults, deformationScale]);

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

      {/* 원래 형태(Original Shape) - 반투명 처리 */}
      <group opacity={0.3} transparent={true}>
        <NodeInstancedMesh nodes={nodes} color={nodeSettings.color} size={nodeSettings.size} />
        {linesGeometry && (
          <lineSegments geometry={linesGeometry}>
            <lineBasicMaterial attach="material" color="gray" opacity={0.4} transparent={true} />
          </lineSegments>
        )}
      </group>

      {/* 노드 텍스트 번호 표시 (원래 위치) */}
      {nodes.map(({ id, x, y, z }, idx) => (
        <NodeLabel key={idx} pos={[x, y, z]} id={id} color="black" size={nodeSettings.size} />
      ))}

      {/* 변형된 형태 (Deformed Shape) - 빨간색으로 오버레이 */}
      {deformedNodes.length > 0 && (
        <group>
          <NodeInstancedMesh nodes={deformedNodes} color="red" size={nodeSettings.size * 1.2} />
          {deformedLinesGeometry && (
            <lineSegments geometry={deformedLinesGeometry}>
              <lineBasicMaterial attach="material" color="red" linewidth={2} />
            </lineSegments>
          )}
        </group>
      )}
    </>
  );
}

export default MyElement3D;
