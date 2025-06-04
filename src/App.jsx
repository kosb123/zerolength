import React, { useState } from 'react';
import './App.css';
import { Canvas } from '@react-three/fiber';
import MyElement3D from './MyElement3D';
import NodeMenu from './Menu/NodeMenu';
import MainMenu from './Menu/MainMenu';
import Members from './Menu/Members';
import SectionMenu from './Menu/SectionMenu';

function App() {
  const [activeMenu, setActiveMenu] = useState('Main');

  function handleMenuClick(menu) {
    if (activeMenu === menu) {
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

      <main id="canvas-container">
        <Canvas>
          <MyElement3D />
        </Canvas>
      </main>
    </div>
  );
}

export default App;
