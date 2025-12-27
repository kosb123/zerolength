// === solve.js (전체 수정 코드) ===

// 1. 필요한 모듈 가져오기
import * as fs from 'fs'; // Node.js 파일 시스템 모듈
import * as THREE from 'three';
import * as math from 'mathjs';

// 2. 데이터 정의 (입력값)
const nodes = [
  { id: 1, x: 0, y: 0, z: 0 },
  { id: 2, x: 0, y: 3, z: 0 },
  { id: 3, x: 3, y: 3, z: 0 },
  { id: 4, x: 6, y: 3, z: 0 },
  { id: 5, x: 9, y: 0, z: 3 }
];

const members = [
  { member_id: 1, n1: 1, n2: 2, sec_id: 1 },
  { member_id: 2, n1: 2, n2: 3, sec_id: 1 },
  { member_id: 3, n1: 3, n2: 4, sec_id: 1 },
  { member_id: 4, n1: 4, n2: 5, sec_id: 1 }
];

const materials = [
  { mat_id: 1, area: 0.01, Iy: 1e-3, Iz: 1e-3, J: 2e-3, E: 200e9, G: 80e9 }
];


const pointloads = [
  { node_id: 3, fx: 0, fy: 0, fz: 240000, mx: 0, my: 0, mz: 0 },
  { node_id: 4, fx: 0, fy: -60000, fz: 0, mx: 0, my: 0, mz: -180000 },
];

const localdistributedloads = [
  { member_id: 1, wy: -40000, wz: 0 }
];


const golbaldistributedloads = [
  { member_id: 1, wy: 0, wz: 0 }
];

// 이 부분은 supports로 대체될 것이므로 제거합니다.
// const displacement = [1, 2, 3, 4, 5, 6, 25, 26, 27, 28, 29, 30];

const supports = [
  { node_id: 1, ux: 0, uy: 0, uz: 0, rx: 0, ry: 0, rz: 0 },
  { node_id: 5, ux: 0, uy: 0, uz: 0, rx: 0, ry: 0, rz: 0 }
];


// 3. 유한요소해석 관련 함수들

// 노드 좌표 차이 계산 함수
function calcNodeDifferences(node1, node2) {
  return {
    dx: node2.x - node1.x,
    dy: node2.y - node1.y,
    dz: node2.z - node1.z
  };
}

// 벡터 길이 계산
function calcLength(dx, dy, dz) {
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

// 요소 로컬 강성 행렬
function createLocalStiffnessMatrix(eal, eiyl, eizl, gjl, el) {
  const K_Axial = eal;
  const K_Torsion = gjl;
  const K_ShearY = 12 * eiyl / (el * el);
  const K_ShearZ = 12 * eizl / (el * el);
  const K_CoupleY = 6 * eiyl / el;
  const K_CoupleZ = 6 * eizl / el;
  const K_Bend4Y = 4 * eiyl;
  const K_Bend2Y = 2 * eiyl;
  const K_Bend4Z = 4 * eizl;
  const K_Bend2Z = 2 * eizl;

  return [
    [K_Axial, 0, 0, 0, 0, 0, -K_Axial, 0, 0, 0, 0, 0],
    [0, K_ShearZ, 0, 0, 0, K_CoupleZ, 0, -K_ShearZ, 0, 0, 0, K_CoupleZ],
    [0, 0, K_ShearY, 0, -K_CoupleY, 0, 0, 0, -K_ShearY, 0, -K_CoupleY, 0],
    [0, 0, 0, K_Torsion, 0, 0, 0, 0, 0, -K_Torsion, 0, 0],
    [0, 0, -K_CoupleY, 0, K_Bend4Y, 0, 0, 0, K_CoupleY, 0, K_Bend2Y, 0],
    [0, K_CoupleZ, 0, 0, 0, K_Bend4Z, 0, -K_CoupleZ, 0, 0, 0, K_Bend2Z],
    [-K_Axial, 0, 0, 0, 0, 0, K_Axial, 0, 0, 0, 0, 0],
    [0, -K_ShearZ, 0, 0, 0, -K_CoupleZ, 0, K_ShearZ, 0, 0, 0, -K_CoupleZ],
    [0, 0, -K_ShearY, 0, K_CoupleY, 0, 0, 0, K_ShearY, 0, K_CoupleY, 0],
    [0, 0, 0, -K_Torsion, 0, 0, 0, 0, 0, K_Torsion, 0, 0],
    [0, 0, -K_CoupleY, 0, K_Bend2Y, 0, 0, 0, K_CoupleY, 0, K_Bend4Y, 0],
    [0, K_CoupleZ, 0, 0, 0, K_Bend2Z, 0, -K_CoupleZ, 0, 0, 0, K_Bend4Z],
  ];
}

// 방향여현 행렬
function getRotationMatrixFromVector(targetVector) {
  const axis_x = new THREE.Vector3(1, 0, 0);
  const axis_y = new THREE.Vector3(0, 1, 0);
  const axis_z = new THREE.Vector3(0, 0, 1);

  const v2 = targetVector.clone().normalize();

  const axis = new THREE.Vector3().crossVectors(axis_x, v2);

  if (axis.length() < 1e-10) {
    if (axis_x.dot(v2) > 0) {
      // 단위 행렬 반환
      return math.identity(12).toArray();
    } else {
      const quaternion = new THREE.Quaternion();
      quaternion.setFromAxisAngle(axis_y, Math.PI);

      const rx = axis_x.clone().applyQuaternion(quaternion);
      const ry = axis_y.clone().applyQuaternion(quaternion);
      const rz = axis_z.clone().applyQuaternion(quaternion);

        return [
    [rx.x, rx.y, rx.z, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [ry.x, ry.y, ry.z, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [rz.x, rz.y, rz.z, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, rx.x, rx.y, rx.z, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, ry.x, ry.y, ry.z, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, rz.x, rz.y, rz.z, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, rx.x, rx.y, rx.z, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, ry.x, ry.y, ry.z, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, rz.x, rz.y, rz.z, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, rx.x, rx.y, rx.z],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, ry.x, ry.y, ry.z],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, rz.x, rz.y, rz.z]
  ];
    }
  }

  axis.normalize();

  const angle = Math.acos(axis_x.dot(v2) / (axis_x.length() * v2.length()));

  const quaternion = new THREE.Quaternion();
  quaternion.setFromAxisAngle(axis, angle);

  const rx = axis_x.clone().applyQuaternion(quaternion);
  const ry = axis_y.clone().applyQuaternion(quaternion);
  const rz = axis_z.clone().applyQuaternion(quaternion);

  return [
    [rx.x, rx.y, rx.z, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [ry.x, ry.y, ry.z, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [rz.x, rz.y, rz.z, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, rx.x, rx.y, rx.z, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, ry.x, ry.y, ry.z, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, rz.x, rz.y, rz.z, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, rx.x, rx.y, rx.z, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, ry.x, ry.y, ry.z, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, rz.x, rz.y, rz.z, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, rx.x, rx.y, rx.z],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, ry.x, ry.y, ry.z],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, rz.x, rz.y, rz.z]
  ];
}

// 요소 전역 강성 행렬 (수정 없음)
function ElemStiff(n) {
  const member = members.find(e => e.member_id === n);
  if (!member) throw new Error(`Element ID ${n} not found!`);
  const node1 = nodes.find(node => node.id === member.n1);
  const node2 = nodes.find(node => node.id === member.n2);
  if (!node1 || !node2) throw new Error("Node(s) not found!");
  const material = materials.find(mat => mat.mat_id === member.sec_id);
  if (!material) throw new Error("Material not found!");

  const { dx, dy, dz } = calcNodeDifferences(node1, node2);
  const L = calcLength(dx, dy, dz);

  const EA_L = material.E * material.area / L;
  const EIy_L = material.E * material.Iy / L;
  const EIz_L = material.E * material.Iz / L;
  const GJ_L = material.G * material.J / L;

  const localStiffnessMatrix = createLocalStiffnessMatrix(EA_L, EIy_L, EIz_L, GJ_L, L);
  const targetVector = new THREE.Vector3(dx, dy, dz);
  const rotationMatrix = getRotationMatrixFromVector(targetVector);

  const tempMatrix = math.multiply(math.transpose(rotationMatrix), localStiffnessMatrix);
  return math.multiply(tempMatrix, rotationMatrix);
}

// 전체 강성 행렬 조립 (수정 없음)
function Stiffness() {
  const maxId = Math.max(...nodes.map(n => n.id));
  const stiff = math.zeros(6 * maxId, 6 * maxId).toArray();

  members.forEach(member => {
    const elemStiff = ElemStiff(member.member_id);
    const n1 = member.n1 - 1;
    const n2 = member.n2 - 1;

    for (let i = 0; i < 12; i++) {
      for (let j = 0; j < 12; j++) {
        const rowIndex = (i < 6) ? n1 * 6 + i : n2 * 6 + (i - 6);
        const colIndex = (j < 6) ? n1 * 6 + j : n2 * 6 + (j - 6);
        stiff[rowIndex][colIndex] += elemStiff[i][j];
      }
    }
  });
  return stiff;
}

// 페널티법 적용 (수정됨)
function PenaltyMethod(stiff, force, supports) { // supports와 nodes를 인수로 받음
    const flatStiff = stiff.flat();
    const maxVal = flatStiff.reduce((max, val) => Math.max(max, Math.abs(val)), 0);
    const cnst = 10000 * maxVal;


    // 각 자유도(ux, uy, uz, rx, ry, rz)의 속성명과 해당 인덱스 오프셋 매핑
    const dofMap = {
        ux: 0,
        uy: 1,
        uz: 2,
        rx: 3,
        ry: 4,
        rz: 5
    };

    supports.forEach(support => {
        const nodeGlobalIndex = (support.node_id - 1) * 6; // 해당 노드의 전역 자유도 시작 인덱스

        // 각 구속 자유도 속성을 순회
        for (const dofKey in dofMap) {

            if (Object.prototype.hasOwnProperty.call(support, dofKey)) {
                if (support[dofKey] === 0) {
                    const idx = nodeGlobalIndex + dofMap[dofKey]; 
                    
                    if (stiff[idx] !== undefined) { 
                        stiff[idx][idx] += cnst;
                        force[idx] += cnst * support[dofKey]; // 또는 force[idx] += cnst * support[dofKey]; 하지만 0이므로 생략 가능
                    }
                }
            }
        }
    });
    return [stiff, force]; // 수정된 강성 행렬을 반환합니다.

}


/**
 * 2차원 배열을 CSV 파일로 저장하는 Node.js용 함수
 * @param {number[][]} matrix - 저장할 2차원 배열
 * @param {string} filename - 저장할 파일 이름
 */
function saveMatrixToCSV(matrix, filename) {
  const csvRows = matrix.map(row => row.join(','));
  const csvString = csvRows.join('\n');

  try {
    fs.writeFileSync(filename, csvString, 'utf8');
    console.log(`✅ ${filename} 파일이 성공적으로 저장되었습니다.`);
  } catch (err) {
    console.error(`❌ 파일 저장 중 오류 발생: ${err}`);
  }
}


function createGlobalForceVector(pointloads, nodes) {
    // 1. nodes 배열에서 최대 node_id 찾기
    // nodes 배열이 비어있을 경우를 대비하여 0을 초기값으로 설정
    const maxNodeId = nodes.length > 0 ? Math.max(...nodes.map(n => n.id)) : 0;

    if (maxNodeId === 0) {
        console.warn("Warning: No nodes found. Global force vector will be empty.");
        return [];
    }

    // 2. 전체 force 배열의 크기 결정 (최대 노드 ID * 각 노드당 자유도 6개)
    const dofsPerNode = 6; // 3D 보 요소: ux, uy, uz, rx, ry, rz
    const totalDofs = maxNodeId * dofsPerNode;

    // 3. force 배열을 0으로 초기화
    const force = new Array(totalDofs).fill(0);

    // 4. pointloads 데이터를 force 배열에 매핑
    pointloads.forEach(load => {
        // node_id가 유효한 범위 내에 있는지 확인 (최대 node_id보다 작거나 같아야 함)
        if (load.node_id < 1 || load.node_id > maxNodeId) {
            console.warn(`Warning: pointload for invalid node_id ${load.node_id}. Skipping.`);
            return; // 해당 pointload 건너뛰기
        }

        const nodeBaseIndex = (load.node_id - 1) * dofsPerNode; // 해당 노드의 전역 자유도 시작 인덱스

        // 각 자유도(fx, fy, fz, mx, my, mz)에 대해 값을 할당
        // 값이 0인 경우에도 명시적으로 0을 할당하여 초기화된 0을 덮어씁니다.
        // 이는 초기화된 0과 명시적으로 0인 하중을 구분하기 위함이지만,
        // 대부분의 경우 `fill(0)`만으로도 충분합니다.
        // 여기서는 `if (load.fx !== 0)`와 같은 조건문을 제거하여 모든 값을 명시적으로 할당합니다.
        
        force[nodeBaseIndex + 0] = load.fx; // ux (0)
        force[nodeBaseIndex + 1] = load.fy; // uy (1)
        force[nodeBaseIndex + 2] = load.fz; // uz (2)
        force[nodeBaseIndex + 3] = load.mx; // rx (3)
        force[nodeBaseIndex + 4] = load.my; // ry (4)
        force[nodeBaseIndex + 5] = load.mz; // rz (5)
    });

    return force;
}

function getLocalDistributedLoads(localdistributedloads, members, nodes, force) {
    // Check if localdistributedloads is empty or not an array
  if (!localdistributedloads || localdistributedloads.length === 0) {
    console.log("No distributed loads to process. Returning early.");
    return force; // 바로 함수 종료
  }

   localdistributedloads.forEach(load => {
    const { member_id, wy, wz } = load;
    
    const memberId = members.find(member => member.member_id === member_id);
      const node1_id = memberId.n1;
      const node2_id = memberId.n2;

      const node1 = nodes.find(node => node.id === node1_id);
      const node2 = nodes.find(node => node.id === node2_id);
      const dx = node2.x - node1.x;
      const dy = node2.y - node1.y;
      const dz = node2.z - node1.z;
      const length = Math.sqrt(dx * dx + dy * dy + dz * dz);

  const targetVector = new THREE.Vector3(dx, dy, dz);
  const rotationMatrix = getRotationMatrixFromVector(targetVector);


      const localLoad = [0, wy*length/2, wz*length/2, 
         0, -wz*length*length/12, wy*length*length/12, 
         0, wy*length/2, wz*length/2,
         0, wz*length*length/12, -wy*length*length/12,];
      const localLoadTransformed = math.multiply(math.inv(rotationMatrix), localLoad);
      
      for (let i = 0; i < 6; i++) {
         force[(node1_id-1) * 6 + i] += localLoadTransformed[i];
         force[(node2_id-1) * 6 + i] += localLoadTransformed[i + 6];

      }
});

 return force;

}



const force1 = createGlobalForceVector(pointloads, nodes);
const force2 = getLocalDistributedLoads(localdistributedloads, members, nodes, force1)




// 4. 메인 실행 로직
console.log("해석을 시작합니다...");

// 1) 전체 강성 행렬 생성 및 저장
const stiff = Stiffness();
saveMatrixToCSV(stiff, 'stiffness_initial.csv');

// 2) 페널티법 적용 및 저장
// supports와 nodes를 인수로 전달
const [stiff_after, force_after] = PenaltyMethod(stiff, force2, supports);
saveMatrixToCSV(stiff_after, 'stiffness_penalty.csv');

// 3) 변위 계산 및 저장
const u = math.multiply(math.inv(stiff_after), force_after);
// 1차원 벡터 u를 CSV로 저장하기 위해 2차원 배열 [u] 형태로 전달
saveMatrixToCSV([u], 'displacement_u.csv');

console.log("해석이 완료되었습니다. 3개의 CSV 파일이 생성되었습니다.");