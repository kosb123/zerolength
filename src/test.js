import * as THREE from 'three';
import * as math from 'mathjs';

const force = [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  240000,  0,  0,  0,  0,  -60000,  0,  0,  0,  -180000,  0,  0,  0,  0,  0,  0,]

const localdistributedloads = [
  { member_id: 1, wy: -4000, wz: 0 },
  { member_id: 2, wy: -4000, wz: 0 }
];

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




function getLocalDistributedLoads(localdistributedloads, members, nodes, force) {
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
getLocalDistributedLoads(localdistributedloads, members, nodes, force);