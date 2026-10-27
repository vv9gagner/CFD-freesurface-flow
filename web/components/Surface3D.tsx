"use client";

import { useEffect, useRef } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

import { colorLut } from "@/lib/colormaps";

type Props = {
  eta: Float32Array; // surface height per cell, nx*ny, x fastest
  nx: number;
  ny: number;
  Lx: number;
  Ly: number;
  Lz: number;
  depth: number;
  range: number; // colour scale covers depth +- range
  exaggeration: number; // vertical scale of eta - depth
  dark: boolean;
};

type Scene3D = {
  renderer: THREE.WebGLRenderer;
  scene: THREE.Scene;
  camera: THREE.PerspectiveCamera;
  controls: OrbitControls;
  surface: THREE.BufferGeometry;
  skirt: THREE.BufferGeometry;
  skirtMaterial: THREE.MeshStandardMaterial;
  boxMaterial: THREE.LineBasicMaterial;
  render: () => void;
  dispose: () => void;
};

/** sRGB lookup table converted to linear RGB, which three.js expects for vertex colours. */
function linearLut(dark: boolean): Float32Array {
  const lut = colorLut("diverging", dark);
  const out = new Float32Array(256 * 3);
  for (let n = 0; n < 256; n++)
    for (let c = 0; c < 3; c++) {
      const v = lut[4 * n + c] / 255;
      out[3 * n + c] = v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
    }
  return out;
}

/** Tank coordinates (x, y, z up) to three.js (X, Y up, Z), centred horizontally; right-handed. */
const toX = (x: number, Lx: number) => x - Lx / 2;
const toZ = (y: number, Ly: number) => Ly / 2 - y;

export default function Surface3D({ eta, nx, ny, Lx, Ly, Lz, depth, range, exaggeration, dark }: Props) {
  const mountRef = useRef<HTMLDivElement>(null);
  const sceneRef = useRef<Scene3D | null>(null);

  // scene, camera, geometry buffers: rebuilt only when the tank or grid changes
  useEffect(() => {
    const mount = mountRef.current;
    if (!mount) return;
    const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
    renderer.setPixelRatio(window.devicePixelRatio || 1);
    renderer.setClearColor(0x000000, 0);
    mount.appendChild(renderer.domElement);
    renderer.domElement.style.display = "block";
    renderer.domElement.style.touchAction = "none";

    const scene = new THREE.Scene();
    scene.add(new THREE.HemisphereLight(0xffffff, 0x444444, 1.4));
    const sun = new THREE.DirectionalLight(0xffffff, 1.6);
    sun.position.set(1, 2.5, 1.5);
    scene.add(sun);

    // free surface: (nx+1) x (ny+1) vertices on the cell corners
    const surface = new THREE.BufferGeometry();
    const nv = (nx + 1) * (ny + 1);
    surface.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3 * nv), 3));
    surface.setAttribute("color", new THREE.BufferAttribute(new Float32Array(3 * nv), 3));
    const index: number[] = [];
    for (let j = 0; j < ny; j++)
      for (let i = 0; i < nx; i++) {
        const a = j * (nx + 1) + i, b = a + 1, c = a + nx + 1, d = c + 1;
        index.push(a, c, b, b, c, d);
      }
    surface.setIndex(index);
    const surfaceMesh = new THREE.Mesh(
      surface,
      new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.5, metalness: 0.05, side: THREE.DoubleSide }),
    );
    scene.add(surfaceMesh);

    // translucent water walls along the four sides, from the bottom up to the surface
    const perimeter = 2 * (nx + 1) + 2 * (ny + 1);
    const skirt = new THREE.BufferGeometry();
    skirt.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3 * 2 * perimeter), 3));
    const sIndex: number[] = [];
    let start = 0;
    for (const len of [nx + 1, nx + 1, ny + 1, ny + 1]) {
      for (let m = 0; m < len - 1; m++) {
        const lo = 2 * (start + m), hi = lo + 1, lo2 = lo + 2, hi2 = lo + 3;
        sIndex.push(lo, lo2, hi, hi, lo2, hi2);
      }
      start += len;
    }
    skirt.setIndex(sIndex);
    const skirtMaterial = new THREE.MeshStandardMaterial({
      transparent: true,
      opacity: 0.28,
      side: THREE.DoubleSide,
      depthWrite: false,
      roughness: 0.8,
    });
    scene.add(new THREE.Mesh(skirt, skirtMaterial));

    // tank outline
    const box = new THREE.EdgesGeometry(new THREE.BoxGeometry(Lx, Lz, Ly));
    const boxMaterial = new THREE.LineBasicMaterial({ transparent: true, opacity: 0.6 });
    const outline = new THREE.LineSegments(box, boxMaterial);
    outline.position.set(0, Lz / 2, 0);
    scene.add(outline);

    const size = Math.max(Lx, Ly, Lz);
    const camera = new THREE.PerspectiveCamera(35, 1, size * 0.01, size * 20);
    camera.position.set(0.85 * size, 1.0 * size, 1.45 * size);
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.target.set(0, depth * 0.6, 0);
    controls.maxPolarAngle = Math.PI * 0.49;
    controls.update();

    const render = () => renderer.render(scene, camera);
    controls.addEventListener("change", render);

    const ro = new ResizeObserver(([entry]) => {
      const w = Math.max(1, Math.floor(entry.contentRect.width));
      const h = Math.min(Math.floor(w * 0.62), Math.floor(window.innerHeight * 0.7));
      renderer.setSize(w, h);
      camera.aspect = w / h;
      camera.updateProjectionMatrix();
      render();
    });
    ro.observe(mount);

    sceneRef.current = {
      renderer,
      scene,
      camera,
      controls,
      surface,
      skirt,
      skirtMaterial,
      boxMaterial,
      render,
      dispose: () => {
        ro.disconnect();
        controls.dispose();
        surface.dispose();
        skirt.dispose();
        box.dispose();
        (surfaceMesh.material as THREE.Material).dispose();
        skirtMaterial.dispose();
        boxMaterial.dispose();
        renderer.dispose();
        renderer.domElement.remove();
      },
    };
    return () => {
      sceneRef.current?.dispose();
      sceneRef.current = null;
    };
  }, [nx, ny, Lx, Ly, Lz, depth]);

  // per-frame update of heights and colours
  useEffect(() => {
    const s = sceneRef.current;
    if (!s || eta.length !== nx * ny) return;
    const dx = Lx / nx, dy = Ly / ny;
    const lut = linearLut(dark);
    const pos = s.surface.getAttribute("position") as THREE.BufferAttribute;
    const col = s.surface.getAttribute("color") as THREE.BufferAttribute;
    const P = pos.array as Float32Array, Cc = col.array as Float32Array;
    const height = (i: number, j: number) => {
      // corner height: mean of the up to four cells that share the corner
      let sum = 0, n = 0;
      for (let jj = j - 1; jj <= j; jj++)
        for (let ii = i - 1; ii <= i; ii++)
          if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
            sum += eta[jj * nx + ii];
            n++;
          }
      return sum / n;
    };
    const Y = (h: number) => depth + (h - depth) * exaggeration;
    const corner = new Float32Array((nx + 1) * (ny + 1));
    for (let j = 0; j <= ny; j++)
      for (let i = 0; i <= nx; i++) {
        const v = j * (nx + 1) + i;
        const h = height(i, j);
        corner[v] = h;
        P[3 * v] = toX(i * dx, Lx);
        P[3 * v + 1] = Y(h);
        P[3 * v + 2] = toZ(j * dy, Ly);
        const t = Math.min(1, Math.max(0, ((h - depth) / (range || 1) + 1) / 2));
        const k = Math.round(t * 255) * 3;
        Cc[3 * v] = lut[k];
        Cc[3 * v + 1] = lut[k + 1];
        Cc[3 * v + 2] = lut[k + 2];
      }
    pos.needsUpdate = true;
    col.needsUpdate = true;
    s.surface.computeVertexNormals();
    s.surface.computeBoundingSphere();

    // side walls follow the boundary of the surface: front (y=0), back (y=Ly), left (x=0), right (x=Lx)
    const SP = (s.skirt.getAttribute("position") as THREE.BufferAttribute).array as Float32Array;
    let w = 0;
    const put = (i: number, j: number) => {
      const X = toX(i * dx, Lx), Z = toZ(j * dy, Ly);
      SP.set([X, 0, Z, X, Y(corner[j * (nx + 1) + i]), Z], 6 * w++);
    };
    for (let i = 0; i <= nx; i++) put(i, 0);
    for (let i = 0; i <= nx; i++) put(i, ny);
    for (let j = 0; j <= ny; j++) put(0, j);
    for (let j = 0; j <= ny; j++) put(nx, j);
    s.skirt.getAttribute("position").needsUpdate = true;
    s.skirt.computeVertexNormals();
    s.skirt.computeBoundingSphere();

    s.skirtMaterial.color.set(dark ? "#3987e5" : "#2a78d6");
    s.boxMaterial.color.set(dark ? "#9a9890" : "#6f6e69");
    s.render();
  }, [eta, nx, ny, Lx, Ly, Lz, depth, range, exaggeration, dark]);

  return (
    <div
      ref={mountRef}
      className="w-full cursor-grab active:cursor-grabbing"
      role="img"
      aria-label="Interactive 3D view of the free surface; drag to rotate, scroll to zoom"
    />
  );
}
