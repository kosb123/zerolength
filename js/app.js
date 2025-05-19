import * as THREE from '../node_modules/three/build/three.module.js';

// 장면
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xffffff);

// 카메라
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.z = 5;


// 렌더러
const renderer = new THREE.WebGLRenderer({antialias: true});
renderer.setSize(window.innerWidth, window.innerHeight);

document.body.appendChild(renderer.domElement);

// 큐브
const geometry = new THREE.BoxGeometry(1, 1, 1);
// const material = new THREE.MeshBasicMaterial({color: 0x44aa88});
const material = new THREE.MeshStandardMaterial({color: 0x44aa88});
const cube = new THREE.Mesh(geometry, material);
cube.position.set(1, 0, 0);
scene.add(cube);


function render(time) {
  time *= 0.0005;  // convert time to seconds
 
  cube.rotation.x = time;
  cube.rotation.y = time;
 
  renderer.render(scene, camera);
 
  requestAnimationFrame(render);
}
requestAnimationFrame(render);