import React, { useState } from 'react';

function SelectionWrapper({ children }) {
  const [selectionBox, setSelectionBox] = useState({
    visible: false,
    startX: 0,
    startY: 0,
    currX: 0,
    currY: 0,
  });

  const handleMouseDown = (e) => {
    // Only trigger on Left Click (button 0)
    if (e.button !== 0) return;

    const container = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - container.left;
    const y = e.clientY - container.top;

    setSelectionBox({
      visible: true,
      startX: x,
      startY: y,
      currX: x,
      currY: y,
    });
  };

  const handleMouseMove = (e) => {
    if (!selectionBox.visible) return;

    const container = e.currentTarget.getBoundingClientRect();
    const x = e.clientX - container.left;
    const y = e.clientY - container.top;

    setSelectionBox(prev => ({
      ...prev,
      currX: x,
      currY: y,
    }));
  };

  const handleMouseUp = () => {
    if (selectionBox.visible) {
      // Logic for selecting elements can go here later
      setSelectionBox(prev => ({ ...prev, visible: false }));
    }
  };

  const getBoxStyle = () => {
    const { startX, startY, currX, currY } = selectionBox;
    const left = Math.min(startX, currX);
    const top = Math.min(startY, currY);
    const width = Math.abs(currX - startX);
    const height = Math.abs(currY - startY);

    return {
      left: left,
      top: top,
      width: width,
      height: height,
      position: 'absolute',
      pointerEvents: 'none', // prevent the box itself from catching pointer events
    };
  };

  return (
    <div
      style={{ position: 'relative', width: '100%', height: '100%' }}
      onPointerDown={handleMouseDown}
      onPointerMove={handleMouseMove}
      onPointerUp={handleMouseUp}
      onPointerLeave={handleMouseUp}
    >
      {children}
      {selectionBox.visible && (
        <div
          className={`selection-box ${selectionBox.currX < selectionBox.startX ? 'crossing' : ''}`}
          style={getBoxStyle()}
        />
      )}
    </div>
  );
}

export default SelectionWrapper;
