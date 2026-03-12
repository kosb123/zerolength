export function makeDraggable(el) {
  let isDragging = false;
  let startX = 0, startY = 0;
  let initialLeft = 0, initialTop = 0;

  el.style.cursor = 'grab';
  el.style.pointerEvents = 'auto';

  const onMouseDown = (e) => {
    // Left click only
    if (e.button !== 0) return;
    
    isDragging = true;
    startX = e.clientX;
    startY = e.clientY;
    
    // Get current bounding rectangle
    const rect = el.getBoundingClientRect();
    const parentRect = el.parentElement 
      ? el.parentElement.getBoundingClientRect() 
      : { left: 0, top: 0 };
    
    // Calculate position relative to its containing block (parent or body)
    initialLeft = rect.left - parentRect.left;
    initialTop = rect.top - parentRect.top;
    
    // Switch from right/bottom to left/top fixed pixel positioning
    el.style.right = 'auto';
    el.style.bottom = 'auto';
    el.style.left = initialLeft + 'px';
    el.style.top = initialTop + 'px';
    
    el.style.cursor = 'grabbing';
    e.preventDefault(); 
    e.stopPropagation();
  };

  const onMouseMove = (e) => {
    if (!isDragging) return;
    const dx = e.clientX - startX;
    const dy = e.clientY - startY;
    
    el.style.left = (initialLeft + dx) + 'px';
    el.style.top = (initialTop + dy) + 'px';
  };

  const onMouseUp = () => {
    if (isDragging) {
      isDragging = false;
      el.style.cursor = 'grab';
    }
  };

  el.addEventListener('mousedown', onMouseDown);
  // Bind move/up to window so dragging continues even if cursor leaves element bounds
  window.addEventListener('mousemove', onMouseMove);
  window.addEventListener('mouseup', onMouseUp);

  return () => {
    el.removeEventListener('mousedown', onMouseDown);
    window.removeEventListener('mousemove', onMouseMove);
    window.removeEventListener('mouseup', onMouseUp);
  };
}
