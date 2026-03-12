import React, { useState, useEffect, useRef } from 'react';
import './App.css';
import { Canvas } from '@react-three/fiber';
import { StatsGl } from '@react-three/drei';
import MyElement3D from './MyElement3D';
import PerformanceMonitor from './PerformanceMonitor';
import NodeMenu from './Menu/NodeMenu';
import MainMenu from './Menu/MainMenu';
import Members from './Menu/Members';
import SectionMenu from './Menu/SectionMenu';
import SettingsMenu from './Menu/SettingsMenu';
import SelectionWrapper from './SelectionWrapper';
import useStore from './store/store';
import { computeDisplacements } from './solver/frameSolver';
import { makeDraggable } from './utils/draggable';

function App() {
  const [activeMenu, setActiveMenu] = useState('Main');

  const setNodes = useStore((state) => state.setNodes);
  const setMembers = useStore((state) => state.setMembers);
  const setMaterials = useStore((state) => state.setMaterials);
  const setConstraints = useStore((state) => state.setConstraints);
  const setForces = useStore((state) => state.setForces);
  const setAnalysisResults = useStore((state) => state.setAnalysisResults);

  const showStatsPanel = useStore((state) => state.showStatsPanel);
  const showEnginePanel = useStore((state) => state.showEnginePanel);

  const statsContainerRef = useRef(null);

  useEffect(() => {
    let cleanup = null;
    if (showStatsPanel && statsContainerRef.current) {
      cleanup = makeDraggable(statsContainerRef.current);
    }
    return () => {
      if (cleanup) cleanup();
    };
  }, [showStatsPanel]);

  useEffect(() => {
    // 3D Box Frame Example (8 Nodes, 12 Members)
    // Y is typically the vertical axis in Three.js
    const nodes = [
      // Base nodes (y = 0)
      { id: 1, x: 0, y: 0, z: 0 },
      { id: 2, x: 5, y: 0, z: 0 },
      { id: 3, x: 5, y: 0, z: 5 },
      { id: 4, x: 0, y: 0, z: 5 },
      // Top nodes (y = 5)
      { id: 5, x: 0, y: 5, z: 0 },
      { id: 6, x: 5, y: 5, z: 0 },
      { id: 7, x: 5, y: 5, z: 5 },
      { id: 8, x: 0, y: 5, z: 5 },
    ];

    const members = [
      // Base beams
      { elem_id: 1, n1: 1, n2: 2, sec_id: 1 },
      { elem_id: 2, n1: 2, n2: 3, sec_id: 1 },
      { elem_id: 3, n1: 3, n2: 4, sec_id: 1 },
      { elem_id: 4, n1: 4, n2: 1, sec_id: 1 },
      // Columns
      { elem_id: 5, n1: 1, n2: 5, sec_id: 1 },
      { elem_id: 6, n1: 2, n2: 6, sec_id: 1 },
      { elem_id: 7, n1: 3, n2: 7, sec_id: 1 },
      { elem_id: 8, n1: 4, n2: 8, sec_id: 1 },
      // Top beams
      { elem_id: 9, n1: 5, n2: 6, sec_id: 1 },
      { elem_id: 10, n1: 6, n2: 7, sec_id: 1 },
      { elem_id: 11, n1: 7, n2: 8, sec_id: 1 },
      { elem_id: 12, n1: 8, n2: 5, sec_id: 1 }
    ];

    const materials = [
      { mat_id: 1, area: 0.1, Iy: 1e-4, Iz: 1e-4, J: 2e-4, E: 200e9, G: 80e9 }
    ];

    // Force vector: 8 nodes * 6 DOFs = 48
    const force = new Array(48).fill(0);
    // Apply lateral force (Fx) and downward gravity load (Fy) on top nodes
    const Fx = 100000;
    const Fy = -500000;
    
    // Node 5 (idx 4, DOFs 24-29)
    force[4 * 6 + 0] = Fx; force[4 * 6 + 1] = Fy;
    // Node 6 (idx 5, DOFs 30-35)
    force[5 * 6 + 0] = Fx; force[5 * 6 + 1] = Fy;
    // Node 7 (idx 6, DOFs 36-41)
    force[6 * 6 + 0] = Fx; force[6 * 6 + 1] = Fy;
    // Node 8 (idx 7, DOFs 42-47)
    force[7 * 6 + 0] = Fx; force[7 * 6 + 1] = Fy;

    // Fixed Supports at Nodes 1, 2, 3, 4 (DOFs 1 through 24)
    const constraints = Array.from({ length: 24 }, (_, i) => i + 1);

    setNodes(nodes);
    setMembers(members);
    setMaterials(materials);
    setForces(force);
    setConstraints(constraints);
  }, []);

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
        {activeMenu === 'Settings' && (
          <SettingsMenu onBack={() => setActiveMenu('Main')} />
        )}
      </nav>

      <main id="canvas-container" style={{ position: 'relative' }}>
        {showStatsPanel && (
          <div
            ref={statsContainerRef}
            style={{
              position: 'absolute',
              top: 0,
              right: 0,
              zIndex: 1000,
              pointerEvents: 'auto',
            }}
          />
        )}
        <SelectionWrapper>
          <Canvas>
            <MyElement3D />
            {showEnginePanel && <PerformanceMonitor />}
            {showStatsPanel && <StatsGl parent={statsContainerRef} />}
          </Canvas>
        </SelectionWrapper>
      </main>
    </div>
  );
}

export default App;
