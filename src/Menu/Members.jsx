import React, { useState, useEffect } from 'react';
import parameter from '../store/store';

function Members({ onBack }) {
  const [elemId, setElemId] = useState(1);
  const [node1, setNode1] = useState(0);
  const [node2, setNode2] = useState(0);
  const [matId, setMatId] = useState(0);

  const elements = parameter(state => state.elements);
  const setElements = parameter(state => state.setElements);

  const [selectedIndex, setSelectedIndex] = useState(null);
  const [hoveredIndex, setHoveredIndex] = useState(null);

  const handleElemIdChange = (e) => {
    const val = e.target.value;
    if (val === '') {
      setElemId('');
      return;
    }
    const num = Number(val);
    if (!isNaN(num) && num >= 1 && Number.isInteger(num)) {
      setElemId(num);
    }
  };

  const handleApply = () => {
    if (elemId === '' || elemId < 1) {
      alert('Element ID는 1 이상의 자연수여야 합니다.');
      return;
    }

    const existingIndex = elements.findIndex(el => el.elem_id === elemId);
    if (existingIndex !== -1) {
      const newArr = [...elements];
      newArr[existingIndex] = { elem_id: elemId, n1: node1, n2: node2, mat_id: matId };
      setElements(newArr);
    } else {
      const newArr = [...elements, { elem_id: elemId, n1: node1, n2: node2, mat_id: matId }];
      setElements(newArr);
    }

    let nextId = elemId + 1;
    const ids = new Set(elements.map(e => e.elem_id));
    while (ids.has(nextId)) {
      nextId += 1;
    }

    setElemId(nextId);
    setNode1(0);
    setNode2(0);
    setMatId(0);
    setSelectedIndex(null);
  };

  const handleDelete = (index) => {
    const newArr = elements.filter((_, idx) => idx !== index);
    setElements(newArr);
    setSelectedIndex(null);
  };

  const handleRowClick = (index) => {
    setSelectedIndex(index);
    const elem = elements[index];
    setElemId(elem.elem_id);
    setNode1(elem.n1);
    setNode2(elem.n2);
    setMatId(elem.mat_id);
  };

  useEffect(() => {
    const elem = elements.find(e => e.elem_id === elemId);
    if (elem) {
      setNode1(elem.n1);
      setNode2(elem.n2);
      setMatId(elem.mat_id);
      setSelectedIndex(elements.indexOf(elem));
    } else {
      setNode1(0);
      setNode2(0);
      setMatId(0);
      setSelectedIndex(null);
    }
  }, [elemId, elements]);

  const sortedElements = [...elements].sort((a, b) => a.elem_id - b.elem_id);

  return (
    <div style={styles.container}>
      <div style={styles.topBar}>
        <button onClick={onBack} style={styles.btn}>← Back</button>
        <button style={styles.btn}>☰</button>
      </div>

      <label style={styles.label}>Element ID (1 이상의 자연수)</label>
      <input
        type="number"
        value={elemId}
        min="1"
        onChange={handleElemIdChange}
        style={styles.input}
      />

      <label style={styles.label}>Node1</label>
      <input
        type="number"
        value={node1}
        onChange={e => setNode1(Number(e.target.value))}
        style={styles.input}
      />

      <label style={styles.label}>Node2</label>
      <input
        type="number"
        value={node2}
        onChange={e => setNode2(Number(e.target.value))}
        style={styles.input}
      />

      <label style={styles.label}>Material ID</label>
      <input
        type="number"
        value={matId}
        onChange={e => setMatId(Number(e.target.value))}
        style={styles.input}
      />

      <div style={styles.buttonRow}>
        <button style={styles.btn}>Help</button>
        <button style={styles.btn} onClick={handleApply}>
          {selectedIndex !== null ? '수정' : '추가'}
        </button>
      </div>

      <div style={styles.savedTable}>
        <strong style={{ color: '#ccc' }}>zustand에 저장된 요소 (elements 배열):</strong>
        <table style={styles.table}>
          <thead>
            <tr>
              <th style={styles.th}>ID</th>
              <th style={styles.th}>Node1</th>
              <th style={styles.th}>Node2</th>
              <th style={styles.th}>MatID</th>
              <th style={styles.th}>삭제</th>
            </tr>
          </thead>
          <tbody>
            {sortedElements.length > 0 ? (
              sortedElements.map(({ elem_id, n1, n2, mat_id }, i) => {
                const isExactSelected = elem_id === elemId;
                return (
                  <tr
                    key={elem_id}
                    style={{
                      backgroundColor: isExactSelected ? '#475569' : 'transparent',
                      cursor: 'pointer',
                    }}
                    onClick={() => handleRowClick(elements.findIndex(e => e.elem_id === elem_id))}
                    onMouseEnter={() => setHoveredIndex(i)}
                    onMouseLeave={() => setHoveredIndex(null)}
                  >
                    <td style={styles.td}>{elem_id}</td>
                    <td style={styles.td}>{n1}</td>
                    <td style={styles.td}>{n2}</td>
                    <td style={styles.td}>{mat_id}</td>
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
                          handleDelete(elements.findIndex(el => el.elem_id === elem_id));
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
                  아직 저장된 요소가 없습니다.
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

export default Members;
