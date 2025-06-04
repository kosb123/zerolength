import React, { useState, useEffect } from 'react';
import parameter from '../store/store';

function SectionMenu({ onBack }) {
  const [secId, setSecId] = useState(1);
  const [name, setName] = useState('');
  const [matId, setMatId] = useState(0);
  const [area, setArea] = useState(0);
  const [Iy, setIy] = useState(0);
  const [Iz, setIz] = useState(0);
  const [J, setJ] = useState(0);

  const sections = parameter(state => state.sections);
  const setSections = parameter(state => state.setSections);

  const [selectedIndex, setSelectedIndex] = useState(null);
  const [hoveredIndex, setHoveredIndex] = useState(null);

  const handleIdChange = (e) => {
    const val = e.target.value;
    if (val === '') {
      setSecId('');
      return;
    }
    const num = Number(val);
    if (!isNaN(num) && num >= 1 && Number.isInteger(num)) {
      setSecId(num);
    }
  };

  const handleApply = () => {
    if (secId === '' || secId < 1) {
      alert('Section ID는 1 이상의 자연수여야 합니다.');
      return;
    }

    const existingIndex = sections.findIndex(s => s.sec_id === secId);
    const newData = { sec_id: secId, name, mat_id: matId, A: area, Iy, Iz, J };
    if (existingIndex !== -1) {
      const newArr = [...sections];
      newArr[existingIndex] = newData;
      setSections(newArr);
    } else {
      const newArr = [...sections, newData];
      setSections(newArr);
    }

    let nextId = secId + 1;
    const ids = new Set(sections.map(s => s.sec_id));
    while (ids.has(nextId)) {
      nextId += 1;
    }

    setSecId(nextId);
    setName('');
    setMatId(0);
    setArea(0);
    setIy(0);
    setIz(0);
    setJ(0);
    setSelectedIndex(null);
  };

  const handleDelete = (index) => {
    const newArr = sections.filter((_, idx) => idx !== index);
    setSections(newArr);
    setSelectedIndex(null);
  };

  const handleRowClick = (index) => {
    setSelectedIndex(index);
    const s = sections[index];
    setSecId(s.sec_id);
    setName(s.name);
    setMatId(s.mat_id);
    setArea(s.A);
    setIy(s.Iy);
    setIz(s.Iz);
    setJ(s.J);
  };

  useEffect(() => {
    const s = sections.find(sec => sec.sec_id === secId);
    if (s) {
      setName(s.name);
      setMatId(s.mat_id);
      setArea(s.A);
      setIy(s.Iy);
      setIz(s.Iz);
      setJ(s.J);
      setSelectedIndex(sections.indexOf(s));
    } else {
      setName('');
      setMatId(0);
      setArea(0);
      setIy(0);
      setIz(0);
      setJ(0);
      setSelectedIndex(null);
    }
  }, [secId, sections]);

  const sortedSections = [...sections].sort((a, b) => a.sec_id - b.sec_id);

  return (
    <div style={styles.container}>
      <div style={styles.topBar}>
        <button onClick={onBack} style={styles.btn}>← Back</button>
        <button style={styles.btn}>☰</button>
      </div>

      <label style={styles.label}>Section ID (1 이상의 자연수)</label>
      <input
        type="number"
        value={secId}
        min="1"
        onChange={handleIdChange}
        style={styles.input}
      />

      <label style={styles.label}>Name</label>
      <input
        type="text"
        value={name}
        onChange={e => setName(e.target.value)}
        style={styles.input}
      />

      <label style={styles.label}>Material ID</label>
      <input
        type="number"
        value={matId}
        onChange={e => setMatId(Number(e.target.value))}
        style={styles.input}
      />

      <label style={styles.label}>A</label>
      <input
        type="number"
        value={area}
        onChange={e => setArea(Number(e.target.value))}
        style={styles.input}
      />

      <label style={styles.label}>Iy</label>
      <input
        type="number"
        value={Iy}
        onChange={e => setIy(Number(e.target.value))}
        style={styles.input}
      />

      <label style={styles.label}>Iz</label>
      <input
        type="number"
        value={Iz}
        onChange={e => setIz(Number(e.target.value))}
        style={styles.input}
      />

      <label style={styles.label}>J</label>
      <input
        type="number"
        value={J}
        onChange={e => setJ(Number(e.target.value))}
        style={styles.input}
      />

      <div style={styles.buttonRow}>
        <button style={styles.btn}>Help</button>
        <button style={styles.btn} onClick={handleApply}>
          {selectedIndex !== null ? '수정' : '추가'}
        </button>
      </div>

      <div style={styles.savedTable}>
        <strong style={{ color: '#ccc' }}>zustand에 저장된 섹션 (sections 배열):</strong>
        <table style={styles.table}>
          <thead>
            <tr>
              <th style={styles.th}>ID</th>
              <th style={styles.th}>Name</th>
              <th style={styles.th}>MatID</th>
              <th style={styles.th}>A</th>
              <th style={styles.th}>Iy</th>
              <th style={styles.th}>Iz</th>
              <th style={styles.th}>J</th>
              <th style={styles.th}>삭제</th>
            </tr>
          </thead>
          <tbody>
            {sortedSections.length > 0 ? (
              sortedSections.map((sec, i) => {
                const isExactSelected = sec.sec_id === secId;
                return (
                  <tr
                    key={sec.sec_id}
                    style={{
                      backgroundColor: isExactSelected ? '#475569' : 'transparent',
                      cursor: 'pointer',
                    }}
                    onClick={() => handleRowClick(sections.findIndex(s => s.sec_id === sec.sec_id))}
                    onMouseEnter={() => setHoveredIndex(i)}
                    onMouseLeave={() => setHoveredIndex(null)}
                  >
                    <td style={styles.td}>{sec.sec_id}</td>
                    <td style={styles.td}>{sec.name}</td>
                    <td style={styles.td}>{sec.mat_id}</td>
                    <td style={styles.td}>{sec.A}</td>
                    <td style={styles.td}>{sec.Iy}</td>
                    <td style={styles.td}>{sec.Iz}</td>
                    <td style={styles.td}>{sec.J}</td>
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
                          handleDelete(sections.findIndex(s => s.sec_id === sec.sec_id));
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
                <td colSpan={8} style={{ textAlign: 'center', padding: 12, color: '#777' }}>
                  아직 저장된 섹션이 없습니다.
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
  input: {
    flexGrow: 1,
    padding: 8,
    borderRadius: 4,
    border: '1px solid #64748b',
    backgroundColor: '#334155',
    color: '#fff',
    marginBottom: 16,
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

export default SectionMenu;
