import React, { useState } from 'react';
import './App.css';
import { Canvas } from '@react-three/fiber';
import MyElement3D from './MyElement3D';
import NodeMenu from './Menu/NodeMenu';
import MainMenu from './Menu/MainMenu';
import MemberForm from './Menu/Members';

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
        {/* {activeMenu === 'Members' && (
          <MemberForm onBack={() => setActiveMenu('Main')} />
        )} */}
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
