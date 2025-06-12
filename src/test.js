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
    for (let i = 0; i < 4; i++) {
      for (let j = 0; j < 3; j++) {
        for (let k = 0; k < 3; k++) {
          R12x12[i * 3 + j][i * 3 + k] = R3x3[j][k];
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
  for (let i = 0; i < 4; i++) {
    for (let j = 0; j < 3; j++) {
      for (let k = 0; k < 3; k++) {
        R12x12[i * 3 + j][i * 3 + k] = R3x3[j][k];
      }
    }
  }
  return R12x12;
}




function getLocalDistributedLoads(localdistributedloads, members, nodes) {
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

      const localLoad = [0, wy*length/2, wz*length/2, 
         0, -wz*length*length/12, wy*length*length/12, 
         0, wy*length/2, wz*length/2,
         0, wz*length*length/12, -wy*length*length/12,];

      console.log(`Member ID: ${member_id}, Local Distributed Load:`, localLoad);
      });
}
getLocalDistributedLoads(localdistributedloads, members, nodes);