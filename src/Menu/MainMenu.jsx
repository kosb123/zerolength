import React, { useState } from 'react';

function MainMenu({ onMenuClick }) {
  const [hoverIndex, setHoverIndex] = useState(null);

  const items = [
    { id: 'Nodes', label: 'Nodes' },
    { id: 'Members', label: 'Members' },
    { id: 'Support', label: 'Support' },
    { id: 'Materials', label: 'Materials' },
    { id: 'Sections', label: 'Sections' },
    { id: 'Variables', label: 'Variables' },
    { id: 'Load', label: 'Load' },
    { id: 'Solve', label: 'Solve' },
  ];

  return (
    <>
      <h2>Menu</h2>
      <ul style={{ padding: 0 }}>
        {items.map((item, index) => (
          <li
            key={item.id}
            onClick={() => onMenuClick(item.id)}
            onMouseEnter={() => setHoverIndex(index)}
            onMouseLeave={() => setHoverIndex(null)}
            style={{
              listStyle: 'none',
              cursor: 'pointer',
              padding: '6px 10px',
              borderRadius: 4,
              backgroundColor: hoverIndex === index ? '#3b82f6' : 'transparent',
              color: hoverIndex === index ? 'white' : 'inherit',
              transition: 'background-color 0.3s ease',
            }}
          >
            {item.label}
          </li>
        ))}
      </ul>
    </>
  );
}

export default MainMenu;
