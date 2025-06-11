/**
 * pointloads 데이터를 기반으로 전체 하중 벡터(force array)를 생성합니다.
 * force 배열의 크기는 최대 노드 ID와 각 노드당 자유도 개수(6)에 따라 결정됩니다.
 *
 * @param {Array<Object>} pointloads - 각 노드에 적용되는 점 하중 정보 (node_id, fx, fy, fz, mx, my, mz)
 * @param {Array<Object>} nodes - 전체 노드 정보 (id, x, y, z)
 * @returns {Array<number>} - 각 자유도에 해당하는 하중 값을 포함하는 1차원 배열
 */
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

// --- 사용 예시 ---

// 기존 데이터 (main.js 또는 solve.js의 상단에 정의된 변수라고 가정)
const nodes = [
    { id: 1, x: 0, y: 0, z: 0 },
    { id: 2, x: 0, y: 3, z: 0 },
    { id: 3, x: 3, y: 3, z: 0 },
    { id: 4, x: 6, y: 3, z: 0 },
    { id: 5, x: 9, y: 0, z: 3 }
];

const pointloads = [
    { node_id: 3, fx: 0, fy: 0, fz: 240000, mx: 0, my: 0, mz: 0 },
    { node_id: 4, fx: 0, fy: -60000, fz: 0, mx: 0, my: 0, mz: -180000 },
];

// 함수 호출
const globalForceVector = createGlobalForceVector(pointloads, nodes);

console.log("생성된 전역 하중 벡터:", globalForceVector);