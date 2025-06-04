import React, { useState } from 'react';
import parameter from '../store/store';

function MemberForm({ onBack }) {
  const members = parameter(state => state.members || []);
  const setMembers = parameter(state => state.setMembers);
  
  // nextMemberId 계산 시 undefined 방지
  const maxId = members.length > 0 ? Math.max(...members.map(m => m.MemberID || 0)) : 0;
  const nextMemberId = maxId + 1;

  // 숫자 입력 상태는 빈 문자열로 초기화
  const [nodeA, setNodeA] = useState('');
  const [nodeB, setNodeB] = useState('');
  const [material, setMaterial] = useState('');
  const [area, setArea] = useState('');
  const [Iy, setIy] = useState('');
  const [Iz, setIz] = useState('');
  const [J, setJ] = useState('');

  // clearInputs 함수는 유지하지만 호출하지 않음
  const clearInputs = () => {
    setNodeA('');
    setNodeB('');
    setMaterial('');
    setArea('');
    setIy('');
    setIz('');
    setJ('');
  };

  const handleApply = () => {
    // 빈값 및 숫자 유효성 체크
    if (
      nodeA === '' || nodeB === '' || material.trim() === '' || area === '' ||
      Iy === '' || Iz === '' || J === ''
    ) {
      alert('모든 항목을 입력해주세요.');
      return;
    }

    // 숫자 변환 후 NaN 체크
    const nA = Number(nodeA);
    const nB = Number(nodeB);
    const a = Number(area);
    const iy = Number(Iy);
    const iz = Number(Iz);
    const j = Number(J);

    if ([nA, nB, a, iy, iz, j].some(num => isNaN(num))) {
      alert('숫자 입력이 올바르지 않습니다.');
      return;
    }

    const newMember = {
      MemberID: nextMemberId,
      N1: nA,
      N2: nB,
      Mat: material.trim(),
      Area: a,
      Iy: iy,
      Iz: iz,
      J: j,
    };

    setMembers([...members, newMember]);
    alert(`Member ${nextMemberId}가 저장되었습니다.`);

    // 입력 유지 위해 clearInputs 호출 안 함
  };

  const handleDelete = () => {
    if (members.length === 0) {
      alert('삭제할 멤버가 없습니다.');
      return;
    }
    setMembers(members.slice(0, -1));
    alert('마지막 멤버가 삭제되었습니다.');
  };

  return (
    <div style={styles.container}>
      <div style={styles.topBar}>
        <button style={styles.backBtn} onClick={onBack}>← Back</button>
      </div>

      <div style={styles.label}>Member ID (자동)</div>
      <input
        type="number"
        value={nextMemberId}
        disabled
        style={{ ...styles.input, backgroundColor: '#555' }}
      />

      <label style={styles.label}>N1 (Node A)</label>
      <input
        type="number"
        value={nodeA}
        onChange={e => setNodeA(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>N2 (Node B)</label>
      <input
        type="number"
        value={nodeB}
        onChange={e => setNodeB(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>Mat# (재료 번호)</label>
      <input
        type="text"
        value={material}
        onChange={e => setMaterial(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>Area (단면적)</label>
      <input
        type="number"
        value={area}
        onChange={e => setArea(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>Iy (단면 2차 모멘트 y방향)</label>
      <input
        type="number"
        value={Iy}
        onChange={e => setIy(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>Iz (단면 2차 모멘트 z방향)</label>
      <input
        type="number"
        value={Iz}
        onChange={e => setIz(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>J (단면 극관성 모멘트)</label>
      <input
        type="number"
        value={J}
        onChange={e => setJ(e.target.value)}
        style={styles.input}
      />

      <div style={styles.buttonRow}>
        <button style={styles.deleteBtn} onClick={handleDelete}>Delete Last</button>
        <button style={styles.applyBtn} onClick={handleApply}>Apply</button>
      </div>

      <div style={{ marginTop: 20 }}>
        <strong>저장된 Members: {members.length} 개</strong>
      </div>
    </div>
  );
}

const styles = {
  container: {
    backgroundColor: '#19232d',
    color: 'white',
    padding: 16,
    fontFamily: 'Arial, sans-serif',
    width: 320,
    borderRadius: 8,
    boxSizing: 'border-box',
  },
  topBar: {
    display: 'flex',
    justifyContent: 'flex-start',
    marginBottom: 12,
  },
  backBtn: {
    backgroundColor: '#444',
    border: 'none',
    borderRadius: 4,
    padding: '6px 12px',
    color: 'white',
    cursor: 'pointer',
  },
  label: {
    fontWeight: 'bold',
    fontSize: 14,
    marginTop: 8,
    marginBottom: 4,
  },
  input: {
    width: '100%',
    padding: 6,
    borderRadius: 4,
    border: 'none',
    marginBottom: 8,
    boxSizing: 'border-box',
  },
  buttonRow: {
    display: 'flex',
    justifyContent: 'space-between',
  },
  deleteBtn: {
    backgroundColor: '#b91c1c',
    border: 'none',
    borderRadius: 4,
    color: 'white',
    padding: '8px 16px',
    cursor: 'pointer',
  },
  applyBtn: {
    backgroundColor: '#2563eb',
    border: 'none',
    borderRadius: 4,
    color: 'white',
    padding: '8px 16px',
    cursor: 'pointer',
  },
};

export default MemberForm;
