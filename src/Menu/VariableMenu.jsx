import React, { useState } from "react";
import parameter from "../store/store";

function VariableMenu({ onBack }) {
  const variables = parameter((state) => state.variables);
  const addVariable = parameter((state) => state.addVariable);

  const [showCreate, setShowCreate] = useState(false);
  const [selectedVar, setSelectedVar] = useState(null);

  const [varName, setVarName] = useState("");
  const [keyInput, setKeyInput] = useState("");
  const [keys, setKeys] = useState(["id"]);
  const [rows, setRows] = useState([]);

  const handleAddRow = () => {
    const nextId = rows.length + 1;
    setRows([...rows, { id: nextId }]);
  };

  const handleRowChange = (rowIndex, key, value) => {
    const newRows = [...rows];
    newRows[rowIndex] = { ...newRows[rowIndex], [key]: value };
    setRows(newRows);
  };

  const handleCreate = () => {
    const trimmed = keyInput
      .split(",")
      .map((s) => s.trim())
      .filter(Boolean);
    const nextKeys = ["id", ...trimmed];
    setKeys(nextKeys);

    // 기본 입력 행 하나 생성
    if (rows.length === 0) {
      setRows([{ id: 1 }]);
    }
  };

  const handleSave = () => {
    addVariable({ name: varName, keys, rows });
    setShowCreate(false);
    setVarName("");
    setKeyInput("");
    setKeys(["id"]);
    setRows([]);
  };

  return (
    <div style={styles.container}>
      <div style={styles.topBar}>
        <button onClick={onBack} style={styles.btn}>
          ← Back
        </button>
        <button style={styles.btn}>☰</button>
      </div>
      <h3>Variables</h3>
      <button
        onClick={() => setShowCreate(true)}
        style={{ ...styles.btn, marginBottom: 10 }}
      >
        만들기
      </button>
      <ul style={{ listStyle: "none", padding: 0 }}>
        {variables.map((v, i) => (
          <li
            key={i}
            onClick={() => setSelectedVar(v)}
            style={{ cursor: "pointer", padding: "4px 0" }}
          >
            {v.name}
          </li>
        ))}
      </ul>

      {showCreate && (
        <div style={styles.modalOverlay}>
          <div style={styles.modal}>
            <h4>새 변수 만들기</h4>
            <label style={styles.label}>이름</label>
            <input
              style={styles.input}
              value={varName}
              onChange={(e) => setVarName(e.target.value)}
            />
            <label style={styles.label}>키 (comma separated)</label>
            <input
              style={styles.input}
              value={keyInput}
              onChange={(e) => setKeyInput(e.target.value)}
            />
            <button
              onClick={handleCreate}
              style={{ ...styles.btn, margin: "8px 0" }}
            >
              테이블 생성
            </button>
            {keys.length > 1 && (
              <>
                <table style={styles.table}>
                  <thead>
                    <tr>
                      {keys.map((k) => (
                        <th key={k} style={styles.th}>
                          {k}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {rows.map((row, rIdx) => (
                      <tr key={rIdx}>
                        {keys.map((k) => (
                          <td key={k} style={styles.td}>
                            <input
                              style={styles.cellInput}
                              value={row[k] || ""}
                              onChange={(e) =>
                                handleRowChange(rIdx, k, e.target.value)
                              }
                            />
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
                <button
                  onClick={handleAddRow}
                  style={{ ...styles.btn, marginTop: 8 }}
                >
                  행 추가
                </button>
                <div style={{ marginTop: 10 }}>
                  <button onClick={handleSave} style={styles.btn}>
                    변수설정하기
                  </button>
                  <button
                    onClick={() => setShowCreate(false)}
                    style={{ ...styles.btn, marginLeft: 8 }}
                  >
                    닫기
                  </button>
                </div>
              </>
            )}
          </div>
        </div>
      )}

      {selectedVar && (
        <div style={styles.modalOverlay}>
          <div style={styles.modal}>
            <h4>{selectedVar.name}</h4>
            <table style={styles.table}>
              <thead>
                <tr>
                  {selectedVar.keys.map((k) => (
                    <th key={k} style={styles.th}>
                      {k}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {selectedVar.rows.map((row, rIdx) => (
                  <tr key={rIdx}>
                    {selectedVar.keys.map((k) => (
                      <td key={k} style={styles.td}>
                        {row[k]}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
            <div style={{ marginTop: 10 }}>
              <button onClick={() => setSelectedVar(null)} style={styles.btn}>
                닫기
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

const baseBtn = {
  backgroundColor: "#64748b",
  border: "none",
  padding: "8px 16px",
  borderRadius: "4px",
  color: "#fff",
  cursor: "pointer",
};

const styles = {
  container: {
    backgroundColor: "#000",
    color: "#fff",
    padding: 20,
    fontFamily: "Arial",
    width: 300,
    height: "100%",
    boxSizing: "border-box",
  },
  topBar: {
    display: "flex",
    justifyContent: "space-between",
    marginBottom: 20,
  },
  btn: baseBtn,
  label: {
    fontSize: 14,
    marginBottom: 6,
    display: "block",
  },
  input: {
    padding: 6,
    borderRadius: 4,
    border: "1px solid #64748b",
    backgroundColor: "#334155",
    color: "#fff",
    width: "100%",
    marginBottom: 8,
    boxSizing: "border-box",
  },
  cellInput: {
    padding: 4,
    border: "none",
    backgroundColor: "transparent",
    color: "#fff",
    width: "100%",
    outline: "none",
    textAlign: "center",
  },
  table: {
    width: "100%",
    borderCollapse: "collapse",
    textAlign: "center",
    marginTop: 8,
  },
  th: {
    border: "1px solid #555",
    padding: "6px 4px",
    backgroundColor: "#1e293b",
    fontSize: 13,
    color: "#ccc",
  },
  td: {
    border: "1px solid #444",
    padding: "4px",
    fontSize: 13,
    color: "#eee",
  },
  modalOverlay: {
    position: "fixed",
    top: 0,
    left: 0,
    width: "100%",
    height: "100%",
    backgroundColor: "rgba(0,0,0,0.6)",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
  },
  modal: {
    backgroundColor: "#000",
    padding: 20,
    borderRadius: 4,
    width: "90%",
    maxWidth: 800,
  },
};

export default VariableMenu;
