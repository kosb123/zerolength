import React, { useState, useEffect } from 'react';
import './App.css';
import { Canvas } from '@react-three/fiber';
import MyElement3D from './MyElement3D';
import NodeMenu from './Menu/NodeMenu';
import MainMenu from './Menu/MainMenu';
import Members from './Menu/Members';
import SectionMenu from './Menu/SectionMenu';
import useStore from './store/store';
import { computeDisplacements } from './solver/frameSolver';

function App() {
  const [activeMenu, setActiveMenu] = useState('Main');

  const setNodes = useStore((state) => state.setNodes);
  const setMembers = useStore((state) => state.setMembers);
  const setMaterials = useStore((state) => state.setMaterials);
  const setConstraints = useStore((state) => state.setConstraints);
  const setForces = useStore((state) => state.setForces);
  const setAnalysisResults = useStore((state) => state.setAnalysisResults);

  useEffect(() => {
    // Initial Data Setup (from Frame.js example)
    const nodes = [
      { id: 1, x: 0, y: 0, z: 0 },
      { id: 2, x: 0, y: 3, z: 0 },
      { id: 3, x: 3, y: 3, z: 0 },
      { id: 4, x: 6, y: 3, z: 0 },
      { id: 5, x: 9, y: 0, z: 3 }
    ];

    const members = [
      { elem_id: 1, n1: 1, n2: 2, sec_id: 1 },
      { elem_id: 2, n1: 2, n2: 3, sec_id: 1 },
      { elem_id: 3, n1: 3, n2: 4, sec_id: 1 },
      { elem_id: 4, n1: 4, n2: 5, sec_id: 1 }
    ];

    const materials = [
      { mat_id: 1, area: 0.01, Iy: 1e-3, Iz: 1e-3, J: 2e-3, E: 200e9, G: 80e9 }
    ];

    const force = new Array(30).fill(0);
    // Node 3 (idx 2): Fz = 240000
    force[2 * 6 + 2] = 240000;
    // Node 4 (idx 3): Fy = -60000, Mz = -180000
    force[3 * 6 + 1] = -60000;
    force[3 * 6 + 5] = -180000;

    // Fixed Supports at Node 1 (1-6) and Node 5 (25-30)
    const constraints = [1, 2, 3, 4, 5, 6, 25, 26, 27, 28, 29, 30];

    setNodes(nodes);
    setMembers(members);
    setMaterials(materials);
    setForces(force);
    setConstraints(constraints);
  }, []);

  // Selection Box State
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

    // Get container offset to handle relative coordinates correctly
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

  // Calculate box style
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
    };
  };

  function handleMenuClick(menu) {
    if (menu === 'Solve') {
      const state = useStore.getState();
      const input = {
        nodes: state.nodes,
        members: state.members,
        materials: state.materials,
        force: state.forces,
        constraints: state.constraints,
      };

      const result = computeDisplacements(input);
      if (result) {
        setAnalysisResults(result);
        console.log("Analysis Result:", result);
        alert("Analysis Completed! Check console for details.");
      } else {
        alert("Analysis Failed!");
      }

    } else if (activeMenu === menu) {
      setActiveMenu('Main');
    } else {
      setActiveMenu(menu);
    }
  }

  return (
    <div id="app-container">
      <nav id="sidebar">
        {activeMenu === 'Main' && (
          <MainMenu onMenuClick={handleMenuClick} />
        )}
        {activeMenu === 'Nodes' && (
          <NodeMenu onBack={() => setActiveMenu('Main')} />
        )}
        {activeMenu === 'Members' && (
          <Members onBack={() => setActiveMenu('Main')} />
        )}
        {activeMenu === 'Sections' && (
          <SectionMenu onBack={() => setActiveMenu('Main')} />
        )}
      </nav>

      <main
        id="canvas-container"
        style={{ position: 'relative' }} // Ensure overlay positioning
        onPointerDown={handleMouseDown}
        onPointerMove={handleMouseMove}
        onPointerUp={handleMouseUp}
        onPointerLeave={handleMouseUp} // Stop dragging if left canvas
      >
        <Canvas>
          <MyElement3D />
        </Canvas>

        {selectionBox.visible && (
          <div
            className={`selection-box ${selectionBox.currX < selectionBox.startX ? 'crossing' : ''}`}
            style={getBoxStyle()}
          />
        )}
      </main>
    </div>
  );
}

export default App;
