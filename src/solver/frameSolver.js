import * as THREE from 'three';
import * as math from 'mathjs';

// Calculate coordinate differences between two nodes
function calcNodeDifferences(node1, node2) {
  return {
    dx: node2.x - node1.x,
    dy: node2.y - node1.y,
    dz: node2.z - node1.z,
  };
}

// Compute vector length
function calcLength(dx, dy, dz) {
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

// Create local stiffness matrix for a frame element
function createLocalStiffnessMatrix(eal, eiyl, eizl, gjl, el) {
  const sep = Array.from({ length: 12 }, () => Array(12).fill(0));

  sep[0][0] = eal;
  sep[0][6] = -eal;
  sep[1][1] = 12 * eizl / (el * el);
  sep[1][5] = 6 * eizl / el;
  sep[1][7] = -sep[1][1];
  sep[1][11] = sep[1][5];
  sep[2][2] = 12 * eiyl / (el * el);
  sep[2][4] = -6 * eiyl / el;
  sep[2][8] = -sep[2][2];
  sep[2][10] = sep[2][4];
  sep[3][3] = gjl;
  sep[3][9] = -gjl;
  sep[4][4] = 4 * eiyl;
  sep[4][8] = 6 * eiyl / el;
  sep[4][10] = 2 * eiyl;
  sep[5][5] = 4 * eizl;
  sep[5][7] = -6 * eizl / el;
  sep[5][11] = 2 * eizl;
  sep[6][6] = eal;
  sep[7][7] = 12 * eizl / (el * el);
  sep[7][11] = -6 * eizl / el;
  sep[8][8] = 12 * eiyl / (el * el);
  sep[8][10] = 6 * eiyl / el;
  sep[9][9] = gjl;
  sep[10][10] = 4 * eiyl;
  sep[11][11] = 4 * eizl;

  for (let i = 0; i < 12; i++) {
    for (let j = i; j < 12; j++) {
      sep[j][i] = sep[i][j];
    }
  }

  return sep;
}

// Rotation matrix from local to global coordinates
function getRotationMatrixFromVector(targetVector) {
  const axisX = new THREE.Vector3(1, 0, 0);
  const axisY = new THREE.Vector3(0, 1, 0);
  const axisZ = new THREE.Vector3(0, 0, 1);

  const v2 = targetVector.clone().normalize();
  const axis = new THREE.Vector3().crossVectors(axisX, v2);

  if (axis.length() < 1e-10) {
    if (axisX.dot(v2) > 0) {
      return math.identity(12).toArray();
    }
    const quaternion = new THREE.Quaternion();
    quaternion.setFromAxisAngle(axisY, Math.PI);

    const rx = axisX.clone().applyQuaternion(quaternion);
    const ry = axisY.clone().applyQuaternion(quaternion);
    const rz = axisZ.clone().applyQuaternion(quaternion);

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
      [0, 0, 0, 0, 0, 0, 0, 0, 0, rz.x, rz.y, rz.z],
    ];
  }

  axis.normalize();
  const angle = Math.acos(axisX.dot(v2) / (axisX.length() * v2.length()));
  const quaternion = new THREE.Quaternion();
  quaternion.setFromAxisAngle(axis, angle);

  const rx = axisX.clone().applyQuaternion(quaternion);
  const ry = axisY.clone().applyQuaternion(quaternion);
  const rz = axisZ.clone().applyQuaternion(quaternion);

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
    [0, 0, 0, 0, 0, 0, 0, 0, 0, rz.x, rz.y, rz.z],
  ];
}

// Element stiffness matrix in global coordinates
function elementStiffness(n, nodes, members, materials) {
  const element = members.find((e) => e.elem_id === n);
  if (!element) throw new Error(`Element ${n} not found`);

  const node1 = nodes.find((node) => node.id === element.n1);
  const node2 = nodes.find((node) => node.id === element.n2);
  if (!node1 || !node2) throw new Error('Node not found');

  const material = materials.find((mat) => mat.mat_id === element.sec_id);
  if (!material) throw new Error('Material not found');

  const { dx, dy, dz } = calcNodeDifferences(node1, node2);
  const L = calcLength(dx, dy, dz);

  const EA_L = (material.E * material.area) / L;
  const EIy_L = (material.E * material.Iy) / L;
  const EIz_L = (material.E * material.Iz) / L;
  const GJ_L = (material.G * material.J) / L;

  const local = createLocalStiffnessMatrix(EA_L, EIy_L, EIz_L, GJ_L, L);
  const targetVector = new THREE.Vector3(dx, dy, dz);
  const rotationMatrix = getRotationMatrixFromVector(targetVector);
  const globalStiffness = math.multiply(math.transpose(rotationMatrix), local);
  return math.multiply(globalStiffness, rotationMatrix);
}

// Assemble global stiffness matrix
function assembleStiffness(nodes, members, materials) {
  const maxId = Math.max(...nodes.map((n) => n.id));
  const size = 6 * maxId;
  const stiff = math.zeros(size, size).toArray();

  members.forEach((elem) => {
    const k = elementStiffness(elem.elem_id, nodes, members, materials);
    const n1 = elem.n1 - 1;
    const n2 = elem.n2 - 1;

    for (let i = 0; i < 12; i++) {
      for (let j = 0; j < 12; j++) {
        const row = i < 6 ? n1 * 6 + i : n2 * 6 + (i - 6);
        const col = j < 6 ? n1 * 6 + j : n2 * 6 + (j - 6);
        stiff[row][col] += k[i][j];
      }
    }
  });

  return stiff;
}

// Apply penalty method boundary conditions
function applyPenalty(stiff, force, constraints) {
  const cnst = 10000 * Math.max(...stiff.flat());
  constraints.forEach((dof) => {
    const idx = dof - 1; // 0-based
    stiff[idx][idx] += cnst;
    force[idx] += cnst * 0; // zero displacement
  });
}

// Main solver function
export function computeDisplacements({ nodes, members, materials, force, constraints }) {
  const stiff = assembleStiffness(nodes, members, materials);
  applyPenalty(stiff, force, constraints);
  const u = math.multiply(math.inv(stiff), force);
  return u;
}

