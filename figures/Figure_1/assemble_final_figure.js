#!/usr/bin/env node

// Assemble the manuscript-only Figure 1 composite from the experimental
// procedures (panel A) and canonical bioinformatics workflow (panel B).

const fs = require("fs");
const path = require("path");

function fail(message) {
  process.stderr.write(`${message}\n`);
  process.exit(1);
}

if (process.argv.includes("--help") || process.argv.includes("-h")) {
  process.stdout.write(
    "Usage: node assemble_final_figure.js [OUTPUT.png] [SHARP_MODULE]\n" +
    "\n" +
    "OUTPUT.png defaults to Workflow/output/Figure_1_FinalFigure.png.\n" +
    "SHARP_MODULE defaults to FIGURE_SHARP_MODULE, WORKFLOW_SHARP_MODULE,\n" +
    "or the installed 'sharp' package.\n"
  );
  process.exit(0);
}

const scriptDirectory = __dirname;
const panelAPath = path.join(scriptDirectory, "Experimental_Procedures.png");
const panelBPath = path.join(
  scriptDirectory,
  "Workflow",
  "output",
  "Analysis_Workflow_generated.png"
);
const outputPath = path.resolve(
  process.argv[2] ||
    path.join(scriptDirectory, "Workflow", "output", "Figure_1_FinalFigure.png")
);
const sharpModule =
  process.argv[3] ||
  process.env.FIGURE_SHARP_MODULE ||
  process.env.WORKFLOW_SHARP_MODULE ||
  "sharp";

for (const inputPath of [panelAPath, panelBPath]) {
  if (!fs.existsSync(inputPath)) fail(`Required panel not found: ${inputPath}`);
  if (path.resolve(inputPath) === outputPath) {
    fail(`Output path must differ from panel input: ${outputPath}`);
  }
}

let sharp;
try {
  sharp = require(sharpModule);
} catch (error) {
  fail(`Could not load sharp from '${sharpModule}': ${error.message}`);
}

const density = 300;
const canvasWidth = 3543;
const outerMargin = 45;
const panelWidth = canvasWidth - 2 * outerMargin;
const labelBandHeight = 120;
const panelGap = 60;
const bottomMargin = 45;

async function preparePanel(inputPath, trimWhitespace) {
  let pipeline = sharp(inputPath, { unlimited: true })
    .rotate()
    .flatten({ background: "#ffffff" });

  if (trimWhitespace) {
    pipeline = pipeline.trim({ background: "#ffffff", threshold: 10 });
  }

  return pipeline
    .resize({
      width: panelWidth,
      fit: "inside",
      withoutEnlargement: false,
      kernel: sharp.kernel.lanczos3,
    })
    .png({ compressionLevel: 9, adaptiveFiltering: true })
    .toBuffer({ resolveWithObject: true });
}

function labelOverlay(canvasHeight, panelBStart) {
  return Buffer.from(
    `<svg width="${canvasWidth}" height="${canvasHeight}" ` +
      `xmlns="http://www.w3.org/2000/svg">` +
      `<style>text{font-family:Arial,Liberation Sans,sans-serif;` +
      `font-size:110px;font-weight:700;fill:#000000}</style>` +
      `<text x="18" y="104">A</text>` +
      `<text x="18" y="${panelBStart + 104}">B</text>` +
      `</svg>`
  );
}

async function main() {
  const panelA = await preparePanel(panelAPath, true);
  const panelB = await preparePanel(panelBPath, false);

  const panelATop = labelBandHeight;
  const panelBStart = panelATop + panelA.info.height + panelGap;
  const panelBTop = panelBStart + labelBandHeight;
  const canvasHeight = panelBTop + panelB.info.height + bottomMargin;

  fs.mkdirSync(path.dirname(outputPath), { recursive: true });

  const result = await sharp({
    create: {
      width: canvasWidth,
      height: canvasHeight,
      channels: 3,
      background: "#ffffff",
    },
  })
    .composite([
      { input: panelA.data, left: outerMargin, top: panelATop },
      { input: panelB.data, left: outerMargin, top: panelBTop },
      { input: labelOverlay(canvasHeight, panelBStart), left: 0, top: 0 },
    ])
    .png({ compressionLevel: 9, adaptiveFiltering: true })
    .withMetadata({ density })
    .toFile(outputPath);

  if (result.width !== canvasWidth || result.height !== canvasHeight) {
    fail(
      `Composite dimensions mismatch: ${result.width}x${result.height}; ` +
        `expected ${canvasWidth}x${canvasHeight}.`
    );
  }

  process.stdout.write(
    `Assembled panel A (${panelA.info.width}x${panelA.info.height}) and ` +
      `panel B (${panelB.info.width}x${panelB.info.height}) into ` +
      `${path.basename(outputPath)} (${result.width}x${result.height}, ` +
      `${density} DPI).\n`
  );
}

main().catch((error) => fail(`Figure assembly failed: ${error.message}`));
