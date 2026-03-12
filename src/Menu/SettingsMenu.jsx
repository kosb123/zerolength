import React from 'react';
import useStore from '../store/store';

function SettingsMenu({ onBack }) {
  const nodeSettings = useStore(state => state.nodeSettings);
  const setNodeSettings = useStore(state => state.setNodeSettings);

  const showStatsPanel = useStore(state => state.showStatsPanel);
  const setShowStatsPanel = useStore(state => state.setShowStatsPanel);

  const showEnginePanel = useStore(state => state.showEnginePanel);
  const setShowEnginePanel = useStore(state => state.setShowEnginePanel);

  const handleSizeChange = (e) => {
    setNodeSettings({ ...nodeSettings, size: parseFloat(e.target.value) });
  };

  const handleColorChange = (e) => {
    setNodeSettings({ ...nodeSettings, color: e.target.value });
  };

  return (
    <div style={{ padding: '20px' }}>
      <h2>Settings</h2>
      
      <div style={{ marginBottom: '15px' }}>
        <label style={{ display: 'block', marginBottom: '5px' }}>
          Node Size: {nodeSettings.size}
        </label>
        <input 
          type="range" 
          min="0.1" 
          max="5" 
          step="0.1" 
          value={nodeSettings.size} 
          onChange={handleSizeChange} 
          style={{ width: '100%' }}
        />
      </div>

      <div style={{ marginBottom: '15px' }}>
        <label style={{ display: 'block', marginBottom: '5px' }}>
          Node Color:
        </label>
        <input 
          type="color" 
          value={nodeSettings.color} 
          onChange={handleColorChange} 
          style={{ width: '100%', height: '40px' }}
        />
      </div>

      <div style={{ marginBottom: '15px' }}>
        <h3 style={{ fontSize: '16px', borderBottom: '1px solid #333', paddingBottom: '5px' }}>Performance</h3>
        <label style={{ display: 'flex', alignItems: 'center', marginBottom: '8px', cursor: 'pointer' }}>
          <input 
            type="checkbox" 
            checked={showStatsPanel} 
            onChange={(e) => setShowStatsPanel(e.target.checked)} 
            style={{ marginRight: '8px' }}
          />
          Show FPS & Resource Stats (Top Right)
        </label>

        <label style={{ display: 'flex', alignItems: 'center', cursor: 'pointer' }}>
          <input 
            type="checkbox" 
            checked={showEnginePanel} 
            onChange={(e) => setShowEnginePanel(e.target.checked)} 
            style={{ marginRight: '8px' }}
          />
          Show Engine Metrics (Bottom Right)
        </label>
      </div>

      <button onClick={onBack} style={{ marginTop: '20px', width: '100%', padding: '8px' }}>
        Back
      </button>
    </div>
  );
}

export default SettingsMenu;
