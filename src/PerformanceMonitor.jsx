import React, { useRef, useEffect } from 'react';
import { useFrame, useThree } from '@react-three/fiber';
import { makeDraggable } from './utils/draggable';

export default function PerformanceMonitor() {
  const { gl } = useThree();
  const divRef = useRef(null);
  const cleanupRef = useRef(null);

  useEffect(() => {
    // React 렌더링에 영향을 주지 않기 위해 순수 DOM 노드를 생성하여 오버레이로 띄웁니다.
    const div = document.createElement('div');
    div.style.position = 'absolute';
    div.style.top = '60px'; // 우측 상단 Stats 패널 바로 아래에 위치
    div.style.right = '0px';
    div.style.backgroundColor = 'rgba(0, 0, 0, 0.85)';
    div.style.color = '#00ffcc';
    div.style.padding = '12px';
    div.style.fontFamily = 'monospace';
    div.style.fontSize = '12px';
    div.style.whiteSpace = 'pre';
    div.style.zIndex = '9999';
    div.style.pointerEvents = 'none';
    div.style.borderBottomLeftRadius = '8px';
    div.style.boxShadow = '0 4px 6px rgba(0,0,0,0.3)';
    div.style.lineHeight = '1.5';
    
    document.body.appendChild(div);
    divRef.current = div;
    
    cleanupRef.current = makeDraggable(div);

    return () => {
      if (cleanupRef.current) cleanupRef.current();
      if (div && document.body.contains(div)) {
        document.body.removeChild(div);
      }
    };
  }, []);

  useFrame(() => {
    if (divRef.current) {
      let memStr = '';
      // Chrome 등에서 지원하는 performance.memory API를 통해 자바스크립트 힙 메모리 확인
      if (performance && performance.memory) {
        const used = (performance.memory.usedJSHeapSize / 1024 / 1024).toFixed(1);
        const total = (performance.memory.totalJSHeapSize / 1024 / 1024).toFixed(1);
        const limit = (performance.memory.jsHeapSizeLimit / 1024 / 1024).toFixed(0);
        memStr = `\nJS Memory  : ${used} / ${total} MB\nLimit      : ${limit} MB`;
      } else {
        memStr = `\nJS Memory  : Not Supported (Use Chrome)`;
      }

      // Three.js 렌더링 파이프라인 지표 (최적화에 핵심적인 Draw Call과 삼각형, 형상정보 확인)
      const drawCalls = gl.info.render.calls;
      const triangles = gl.info.render.triangles;
      const geometries = gl.info.memory.geometries;

      divRef.current.innerText = 
        `[ Engine Info ]` + 
        `\nDraw Calls : ${drawCalls}` + 
        `\nTriangles  : ${triangles}` + 
        `\nGeometries : ${geometries}` +
        `\n-----------------------` +
        memStr;
    }
  });

  return null; // 시각적 컴포넌트는 반환하지 않고 로직만 실행
}
