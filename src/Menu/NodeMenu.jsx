// NodeMenu.jsx
import React, { useState, useEffect } from 'react';
import parameter from '../store/store';

function NodeMenu({ onBack }) {
  const [id, setId] = useState(1);
  const [xPos, setXPos] = useState(0);
  const [yPos, setYPos] = useState(0);
  const [zPos, setZPos] = useState(0);

  const nodes = parameter(state => state.nodes);
  const setNodes = parameter(state => state.setNodes);

  const [selectedIndex, setSelectedIndex] = useState(null);
  const [hoveredIndex, setHoveredIndex] = useState(null);

  const handleIdChange = (e) => {
    const val = e.target.value;
    if (val === '') {
      setId('');
      return;
    }
    const num = Number(val);
    if (!isNaN(num) && num >= 1 && Number.isInteger(num)) {
      setId(num);
    }
  };

  const handleApply = () => {
  if (id === '' || id < 1) {
    alert('ID는 1 이상의 자연수여야 합니다.');
    return;
  }

  const existingIndex = nodes.findIndex(node => node.id === id);
  if (existingIndex !== -1) {
    // 수정
    const newArr = [...nodes];
    newArr[existingIndex] = { id, x: xPos, y: yPos, z: zPos };
    setNodes(newArr);
  } else {
    // 추가
    const newArr = [...nodes, { id, x: xPos, y: yPos, z: zPos }];
    setNodes(newArr);
  }

  // 다음 ID 계산: 중복이 없도록 계속 +1
  let nextId = id + 1;
  const ids = new Set(nodes.map(n => n.id));
  while (ids.has(nextId)) {
    nextId += 1;
  }

  setId(nextId);
  setXPos(0);
  setYPos(0);
  setZPos(0);
  setSelectedIndex(null);
};

  const handleDelete = (index) => {
    const newArr = nodes.filter((_, idx) => idx !== index);
    setNodes(newArr);
    setSelectedIndex(null);
  };

  const handleRowClick = (index) => {
    setSelectedIndex(index);
    const node = nodes[index];
    setId(node.id);
    setXPos(node.x);
    setYPos(node.y);
    setZPos(node.z);
  };

  useEffect(() => {
    const node = nodes.find(n => n.id === id);
    if (node) {
      setXPos(node.x);
      setYPos(node.y);
      setZPos(node.z);
      setSelectedIndex(nodes.indexOf(node));
    } else {
      setXPos(0);
      setYPos(0);
      setZPos(0);
      setSelectedIndex(null);
    }
  }, [id, nodes]);

  // 정렬된 노드 배열
  const sortedNodes = [...nodes].sort((a, b) => a.id - b.id);

  return (
    <div style={styles.container}>
      <div style={styles.topBar}>
        <button onClick={onBack} style={styles.btn}>← Back</button>
        <button style={styles.btn}>☰</button>
      </div>

      <label style={styles.label}>Node ID (1 이상의 자연수)</label>
      <input
        type="number"
        value={id}
        min="1"
        onChange={handleIdChange}
        style={styles.input}
      />

      {['X', 'Y', 'Z'].map((axis, i) => {
        const setters = [setXPos, setYPos, setZPos];
        const values = [xPos, yPos, zPos];
        return (
          <div key={axis}>
            <label style={styles.label}>{axis} Position</label>
            <div style={styles.inputRow}>
              <input
                type="number"
                value={values[i]}
                onChange={e => setters[i](Number(e.target.value))}
                style={styles.input}
              />
              <span style={styles.unit}>m</span>
            </div>
          </div>
        );
      })}

      <div style={styles.buttonRow}>
        <button style={styles.btn}>Help</button>
        <button style={styles.btn} onClick={handleApply}>
          {selectedIndex !== null ? '수정' : '추가'}
        </button>
      </div>

      <div style={styles.savedTable}>
        <strong style={{ color: '#ccc' }}>zustand에 저장된 좌표 (nodes 배열):</strong>
        <table style={styles.table}>
          <thead>
            <tr>
              <th style={styles.th}>ID</th>
              <th style={styles.th}>X</th>
              <th style={styles.th}>Y</th>
              <th style={styles.th}>Z</th>
              <th style={styles.th}>삭제</th>
            </tr>
          </thead>
          <tbody>
            {sortedNodes.length > 0 ? (
              sortedNodes.map(({ id: nodeId, x, y, z }, i) => {
                const isExactSelected = nodeId === id;
                return (
                  <tr
                    key={nodeId}
                    style={{
                      backgroundColor: isExactSelected ? '#475569' : 'transparent',
                      cursor: 'pointer',
                    }}
                    onClick={() => handleRowClick(nodes.findIndex(n => n.id === nodeId))}
                    onMouseEnter={() => setHoveredIndex(i)}
                    onMouseLeave={() => setHoveredIndex(null)}
                  >
                    <td style={styles.td}>{nodeId}</td>
                    <td style={styles.td}>{x}</td>
                    <td style={styles.td}>{y}</td>
                    <td style={styles.td}>{z}</td>
                    <td style={styles.td}>
                      <button
                        style={{
                          ...styles.btn,
                          backgroundColor: '#dc2626',
                          padding: '4px 8px',
                          fontSize: 12,
                          visibility: hoveredIndex === i ? 'visible' : 'hidden',
                        }}
                        onClick={e => {
                          e.stopPropagation();
                          handleDelete(nodes.findIndex(n => n.id === nodeId));
                        }}
                      >
                        삭제
                      </button>
                    </td>
                  </tr>
                );
              })
            ) : (
              <tr>
                <td colSpan={5} style={{ textAlign: 'center', padding: 12, color: '#777' }}>
                  아직 저장된 노드가 없습니다.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
}

const baseBtn = {
  backgroundColor: '#64748b',
  border: 'none',
  padding: '8px 16px',
  borderRadius: '4px',
  color: '#fff',
  cursor: 'pointer',
};

const styles = {
  container: {
    backgroundColor: '#000',
    color: '#fff',
    padding: 20,
    fontFamily: 'Arial',
    width: 300,
    height: '100%',
    boxSizing: 'border-box',
  },
  topBar: {
    display: 'flex',
    justifyContent: 'space-between',
    marginBottom: 20,
  },
  btn: baseBtn,
  label: {
    fontSize: 14,
    marginBottom: 6,
    display: 'block',
  },
  inputRow: {
    display: 'flex',
    alignItems: 'center',
    marginBottom: 16,
  },
  input: {
    flexGrow: 1,
    padding: 8,
    borderRadius: 4,
    border: '1px solid #64748b',
    backgroundColor: '#334155',
    color: '#fff',
  },
  unit: {
    marginLeft: 8,
    fontSize: 14,
    color: '#94a3b8',
  },
  buttonRow: {
    display: 'flex',
    justifyContent: 'flex-end',
    gap: 12,
    marginBottom: 20,
  },
  savedTable: {
    marginTop: 10,
  },
  table: {
    width: '100%',
    borderCollapse: 'collapse',
    textAlign: 'center',
    marginTop: 8,
  },
  th: {
    borderBottom: '2px solid #555',
    padding: '8px 4px',
    backgroundColor: '#1e293b',
    fontSize: 13,
    color: '#ccc',
  },
  td: {
    borderBottom: '1px solid #333',
    padding: '6px 4px',
    fontSize: 13,
    color: '#eee',
  },
};

export default NodeMenu;
