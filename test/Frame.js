import * as THREE from '../node_modules/three/build/three.module.js';
import * as math from 'mathjs';


const nodes = [
  { id: 1, x: 0,  y: 0,  z: 0 },
  { id: 2, x: 0,  y: 3,  z: 0 },
  { id: 3, x: 3,  y: 3,  z: 0 },
  { id: 4, x: 6,  y: 3,  z: 0 },
  { id: 5, x: 9,  y: 0,  z: 3 }
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

// const force = [
//   { node_id: 1, fx: 0, fy: 0, fz: 0, mx: 0, my: 0, mz: 0 },
//   { node_id: 2, fx: 0, fy: 0, fz: 0, mx: 0, my: 0, mz: 0 },
//   { node_id: 3, fx: 0, fy: 0, fz: 240000, mx: 0, my: 0, mz: 0 },
//   { node_id: 4, fx: 0, fy: -60000, fz: 0, mx: 0, my: 0, mz: -180000 },
//   { node_id: 5, fx: 0, fy: 0, fz: 0, mx: 0, my: 0, mz: 0 }
// ];

const force = [
  0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0,
  0, 0, 240000, 0, 0, 0,
  0, -60000, 0, 0, 0, -180000,
  0, 0, 0, 0, 0, 0
]


// const displacement = [
//   { node_id: 1, fx: 0, fy: 0, fz: 0, mx: 0, my: 0, mz: 0 },
//   { node_id: 2, fx: null, fy: null, fz: null, mx: null, my: null, mz: null },
//   { node_id: 3, fx: null, fy: null, fz: null, mx: null, my: null, mz: null },
//   { node_id: 4, fx: null, fy: null, fz: null, mx: null, my: null, mz: null },
//   { node_id: 5, fx: 0, fy: 0, fz: 0, mx: 0, my: 0, mz: 0 }
// ];

const displacement = [1, 2, 3, 4, 5, 6, 25, 26, 27, 28, 29, 30];


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
  return Math.sqrt(dx*dx + dy*dy + dz*dz);
}

// 요소 로컬 강성 행렬
// eal: 축방향 강성, eiyl: y축 굽힘 강성, eizl: z축 굽힘 강성, gjl: 전단 강성, el: 길이
function createLocalStiffnessMatrix (eal, eiyl, eizl, gjl, el) {
  let sep = new Array(12);
  for (let i = 0; i < 12; i++) {
    sep[i] = new Array(12).fill(0);  // 2차원 배열 초기화
  }

  sep[0][0] = eal;
  sep[0][6] = -eal;
       sep[1][1] = 12 * eizl / (el*el)
       sep[1][5] = 6 * eizl / el
       sep[1][7] = -sep[1][1]
       sep[1][11] = sep[1][5]
     sep[2][2] = 12 * eiyl / (el*el)
     sep[2][4] = -6 * eiyl / el
     sep[2][8] = -sep[2][2]
     sep[2][10] = sep[2][4]
       sep[3][3] = gjl
       sep[3][9] = -gjl
     sep[4][4] = 4 * eiyl
     sep[4][8] = 6 * eiyl / el
     sep[4][10] = 2 * eiyl
       sep[5][5] = 4 * eizl
       sep[5][7] = -6 * eizl / el
       sep[5][11] = 2 * eizl
     sep[6][6] = eal
       sep[7][7] = 12 * eizl / (el*el)
       sep[7][11] = -6 * eizl / el
     sep[8][8] = 12 * eiyl / (el*el)
     sep[8][10] = 6 * eiyl / el
       sep[9][9] = gjl
     sep[10][10] = 4 * eiyl
       sep[11][11] = 4 * eizl
  for (let i = 0; i < 12; i++) {
    for (let j = i; j < 12; j++) {
      sep[j][i] = sep[i][j];
    }
  }
     
     return sep;
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

// 요소 전역 강성 행렬
function ElemStiff(n) {

    // 요소 정보 저장
  const element = members.find(e => e.elem_id === n);
  if (!element) {
    console.error(`Element ID ${n} not found!`);
    return null;
  }

  const node1 = nodes.find(node => node.id === element.n1);
  const node2 = nodes.find(node => node.id === element.n2);
  if (!node1 || !node2) {
    console.error("Node(s) not found!");
    return null;
  }

  const material = materials.find(mat => mat.mat_id === element.sec_id);
  if (!material) {
    console.error("Material not found!");
    return null;
  }


  // 로컬 강성행렬 구성
  const { dx, dy, dz } = calcNodeDifferences(node1, node2);
  const L = calcLength(dx, dy, dz);

  const EA_L = material.E * material.area / L;
  const EIy_L = material.E * material.Iy / L;
  const EIz_L = material.E * material.Iz / L;
  const GJ_L = material.G * material.J / L;


  const localStiffnessMatrix = createLocalStiffnessMatrix (EA_L, EIy_L, EIz_L, GJ_L, L);
  // console.table(localStiffnessMatrix);


  // 방향여현 행렬 계산
  const targetVector = new THREE.Vector3(node2.x - node1.x, node2.y - node1.y, node2.z - node1.z);
  const rotationMatrix = getRotationMatrixFromVector(targetVector);
  

  // 전역 강성 행렬 계산
  const globalStiffnessMatrix = math.multiply(math.transpose(rotationMatrix), localStiffnessMatrix);
  const finalStiffnessMatrix = math.multiply(globalStiffnessMatrix, rotationMatrix);

  // console.table(rotationMatrix)
  // console.table(localStiffnessMatrix)
  // console.table(finalStiffnessMatrix);

  return finalStiffnessMatrix;
    // localStiffnessMatrix,
    // rotationMatrix,
}

function Stiffness() {

  // 전체 강성행렬 배열 생성
  const maxIdNode = nodes.reduce((maxNode, currentNode) => {
    return currentNode.id > maxNode.id ? currentNode : maxNode;
  });
  const stiff = math.zeros(6 * maxIdNode.id, 6 * maxIdNode.id).toArray();
  // console.table(stiff);

  // 각 요소의 강성 행렬을 전체 강성 행렬에 추가
  members.forEach(element => {
    const elemStiff = ElemStiff(element.elem_id);
      const n1 = element.n1 - 1; // 0-based index
      const n2 = element.n2 - 1; // 0-based index

      for (let i = 0; i < 12; i++) {
        for (let j = 0; j < 12; j++) {
          // 행, 열 인덱스 계산
          const rowIndex = (i < 6) ? n1 * 6 + i : n2 * 6 + (i - 6);
          const colIndex = (j < 6) ? n1 * 6 + j : n2 * 6 + (j - 6);

          // 기존 값에 sourceMatrix[i][j] 더하기
          stiff[rowIndex][colIndex] += elemStiff[i][j];
  }
}
    }
  );
  return stiff;
}

function PenaltyMethod(stiff, force) {
  const cnst = 10000 * Math.max(...stiff.flat());
  
  for (let i = 0; i < displacement.length; i++) {
    const num = displacement[i] - 1; // 0-based index
    stiff[num][num] += cnst;
    force[num] += cnst * 0; // 변위가 0인 경우
  }
}

function downloadMatrixCSV(matrix, filename = 'stiff.csv') {
  // matrix: 2차원 숫자 배열 (예: stiff.toArray())

  const csvRows = matrix.map(row => row.join(',')); // 각 행을 ','로 이어서 문자열로 만듦
  const csvString = csvRows.join('\n');             // 행들을 줄바꿈 문자로 연결

  const blob = new Blob([csvString], { type: 'text/csv' });
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  link.click();
}

const stiff = Stiffness();
console.table(stiff);


PenaltyMethod(stiff, force);
console.table(stiff);

const u = math.multiply(math.inv(stiff), force);
console.table(u);

// 다운로드 호출
// downloadMatrixCSV(stiff)