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
  { elem_id: 1, n1: 1, n2: 2, sec_id: 1 },
  { elem_id: 2, n1: 2, n2: 3, sec_id: 1 },
  { elem_id: 3, n1: 3, n2: 4, sec_id: 1 },
  { elem_id: 4, n1: 4, n2: 5, sec_id: 1 }
];

const materials = [
  { mat_id: 1, area: 0.01, Iy: 1e-3, Iz: 1e-3, J: 2e-3, E: 200e9, G: 80e9 }
];

const force = [
  0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0,
  0, 0, 240000, 0, 0, 0,
  0, -60000, 0, 0, 0, -180000,
  0, 0, 0, 0, 0, 0
];

const displacement = [1, 2, 3, 4, 5, 6, 25, 26, 27, 28, 29, 30];

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
  const v2 = targetVector.clone().normalize();
  const axis = new THREE.Vector3().crossVectors(axis_x, v2);

  if (axis.length() < 1e-10) {
    if (axis_x.dot(v2) > 0) return math.identity(12).toArray();
    
    const quaternion = new THREE.Quaternion().setFromAxisAngle(axis_y, Math.PI);
    const r = new THREE.Matrix4().makeRotationFromQuaternion(quaternion);
    const R3x3 = [
        [r.elements[0], r.elements[1], r.elements[2]],
        [r.elements[4], r.elements[5], r.elements[6]],
        [r.elements[8], r.elements[9], r.elements[10]]
    ];
    const R12x12 = math.zeros(12, 12).toArray();
    for(let i=0; i<4; i++) {
        for(let j=0; j<3; j++) {
            for(let k=0; k<3; k++) {
                R12x12[i*3+j][i*3+k] = R3x3[j][k];
            }
        }
    }
    return R12x12;
  }

  const angle = Math.acos(axis_x.dot(v2));
  const quaternion = new THREE.Quaternion().setFromAxisAngle(axis.normalize(), angle);
  const r = new THREE.Matrix4().makeRotationFromQuaternion(quaternion);
    const R3x3 = [
        [r.elements[0], r.elements[1], r.elements[2]],
        [r.elements[4], r.elements[5], r.elements[6]],
        [r.elements[8], r.elements[9], r.elements[10]]
    ];
    const R12x12 = math.zeros(12, 12).toArray();
    for(let i=0; i<4; i++) {
        for(let j=0; j<3; j++) {
            for(let k=0; k<3; k++) {
                R12x12[i*3+j][i*3+k] = R3x3[j][k];
            }
        }
    }
    return R12x12;
}

// 요소 전역 강성 행렬
function ElemStiff(n) {
  const element = members.find(e => e.elem_id === n);
  if (!element) throw new Error(`Element ID ${n} not found!`);
  const node1 = nodes.find(node => node.id === element.n1);
  const node2 = nodes.find(node => node.id === element.n2);
  if (!node1 || !node2) throw new Error("Node(s) not found!");
  const material = materials.find(mat => mat.mat_id === element.sec_id);
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

// 전체 강성 행렬 조립
function Stiffness() {
  const maxId = Math.max(...nodes.map(n => n.id));
  const stiff = math.zeros(6 * maxId, 6 * maxId).toArray();

  members.forEach(element => {
    const elemStiff = ElemStiff(element.elem_id);
    const n1 = element.n1 - 1;
    const n2 = element.n2 - 1;

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

// 페널티법 적용
function PenaltyMethod(stiff, force) {
  const flatStiff = stiff.flat();
  const maxVal = flatStiff.reduce((max, val) => Math.max(max, Math.abs(val)), 0);
  const cnst = 10000 * maxVal;

  displacement.forEach(dof => {
    const idx = dof - 1;
    if (stiff[idx]) {
      stiff[idx][idx] += cnst;
      // 변위가 0인 경우에 대한 하중은 0을 더하므로 생략 가능
      // force[idx] += cnst * 0;
    }
  });
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


// 4. 메인 실행 로직
console.log("해석을 시작합니다...");

// 1) 전체 강성 행렬 생성 및 저장
const stiff = Stiffness();
saveMatrixToCSV(stiff, 'stiffness_initial.csv');

// 2) 페널티법 적용 및 저장
PenaltyMethod(stiff, force);
saveMatrixToCSV(stiff, 'stiffness_penalty.csv');

// 3) 변위 계산 및 저장
const u = math.multiply(math.inv(stiff), force);
// 1차원 벡터 u를 CSV로 저장하기 위해 2차원 배열 [u] 형태로 전달
saveMatrixToCSV([u], 'displacement_u.csv');

console.log("해석이 완료되었습니다. 3개의 CSV 파일이 생성되었습니다.");