#!/usr/bin/env node

// Rasterize the canonical SVG without reimplementing any visual content.

const fs = require("fs");
const path = require("path");

function fail(message) {
  process.stderr.write(`${message}\n`);
  process.exit(1);
}

if (process.argv.length !== 8) {
  fail("Usage: node rasterize_workflow_svg.js INPUT.svg OUTPUT.png WIDTH HEIGHT DPI SHARP_MODULE");
}

const [, , inputPath, outputPath, widthText, heightText, dpiText, modulePath] = process.argv;
const width = Number.parseInt(widthText, 10);
const height = Number.parseInt(heightText, 10);
const dpi = Number.parseInt(dpiText, 10);

if (!fs.existsSync(inputPath)) fail(`Canonical SVG not found: ${inputPath}`);
if (![width, height, dpi].every(Number.isInteger)) fail("WIDTH, HEIGHT, and DPI must be integers.");
if (width < 1 || height < 1 || dpi < 72) fail("Invalid raster dimensions or density.");

let sharp;
try {
  sharp = require(modulePath);
} catch (error) {
  fail(`Could not load sharp from '${modulePath}': ${error.message}`);
}

fs.mkdirSync(path.dirname(outputPath), { recursive: true });

sharp(inputPath, { density: dpi, unlimited: true })
  .resize(width, height, { fit: "fill", kernel: sharp.kernel.lanczos3 })
  .png({ compressionLevel: 9, adaptiveFiltering: true })
  .withMetadata({ density: dpi })
  .toFile(outputPath)
  .then((info) => {
    if (info.width !== width || info.height !== height) {
      fail(`Raster dimensions mismatch: ${info.width}x${info.height}`);
    }
    process.stdout.write(`Rasterized ${path.basename(inputPath)} -> ${path.basename(outputPath)} (${width}x${height}, ${dpi} DPI)\n`);
  })
  .catch((error) => fail(`SVG rasterization failed: ${error.message}`));
