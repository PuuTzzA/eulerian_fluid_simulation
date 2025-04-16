// import { Fluid } from "./logic.js";

let drawVelX = 0;
let drawVelY = 0;
let drawRed = 0;

let mode = "red";
const modeInidator = document.getElementById("mode");

window.addEventListener("pointerdown", e => {
    if (e.target.id != "canvas") { return };
    const x = Math.floor(e.offsetX / fluid.gridSize);
    const offsetTop = window.innerHeight - fluid.gridSize * fluid.resolutionY;
    const y = Math.floor((e.offsetY - offsetTop) / fluid.gridSize);

    switch (mode) {
        case "velocity":
            fluid.velocitiesX[y * (fluid.resolutionX + 1) + x] = drawVelX;
            fluid.velocitiesX[y * (fluid.resolutionX + 1) + x + 1] = drawVelX;
            fluid.velocitiesY[y * (fluid.resolutionX) + x] = drawVelY;
            fluid.velocitiesY[(y + 1) * (fluid.resolutionX) + x] = drawVelY;
            break;
        case "velocityX":
            fluid.velocitiesX[y * (fluid.resolutionX + 1) + x] = drawVelX;
            fluid.velocitiesX[y * (fluid.resolutionX + 1) + x + 1] = drawVelX;
            break;
        case "velocityY":
            fluid.velocitiesY[y * (fluid.resolutionX) + x] = drawVelY;
            fluid.velocitiesY[(y + 1) * (fluid.resolutionX) + x] = drawVelY;
            break;
        case "red":
            fluid.colorRed[y * fluid.resolutionX + x] = drawRed;
            break;
    }
})

addEventListener("keydown", (event) => {
    if (event.key == "v") {
        mode = "velocity";
    } else if (event.key == "x") {
        mode = "velocityX";
    } else if (event.key == "y") {
        mode = "velocityY";
    } else if (event.key == "r") {
        mode = "red";
    }
    modeInidator.innerHTML = mode;
});

document.addEventListener("mousemove", (e) => {
    //if (e.target.id != "canvas") { return };
    const x = e.offsetX / fluid.gridSize - 0.5;
    const offsetTop = window.innerHeight - fluid.gridSize * fluid.resolutionY;
    const y = (e.offsetY - offsetTop) / fluid.gridSize - 0.5;

    /*     console.log(x, y, "vel", Fluid.bilinearInterpolation(x, y, fluid.colorRed, fluid.resolutionX, fluid.resolutionY));
        if (Fluid.bilinearInterpolation(x, y, fluid.colorRed, fluid.resolutionX, fluid.resolutionY) !=
            Fluid.bilinearInterpolation2(x, y, fluid.colorRed, fluid.resolutionX, fluid.resolutionY)) {
            console.log("bilinear interpolation error");
        } */
    //console.log(x, y, "vel", fluid.getVelocityAtPoint(x, y));
})

window.addEventListener("pointerup", e => {
})

window.addEventListener("resize", e => {
    fluid.initGrid(fluid.resolutionX, fluid.resolutionY);
})

const fps = document.getElementById("fps");
let newFps = 0;
let accumulatedFps = [];
let running = true;
let delta = 0.01;
let previous;

let resolutionX = 1;
let resolutionY = 1;
let density = 0.1;

let fluid = new Fluid(1, 1);

const FIX_DELTA = 0.001;

/* fluid.setResolutionX(5);
fluid.setResolutionY(5);
fluid.update(FIX_DELTA); */
fluid.draw();

function toggleDetailedOptions() {
    document.getElementById("detailed-options").classList.toggle("detailed-options");
    const current = document.getElementById("show-detailed-options").innerHTML;
    document.getElementById("show-detailed-options").innerHTML = current == "arrow_drop_up" ? "arrow_drop_down" : "arrow_drop_up";
}

function toggleDropdown() {
    document.getElementById("quick-select").classList.toggle("dropdown-show");
    const current = document.getElementById("dropdown-button-arrow").innerHTML;
    document.getElementById("dropdown-button-arrow").innerHTML = current == "arrow_drop_up" ? "arrow_drop_down" : "arrow_drop_up";
}

const presets = {
    "Water": [5, 50, 50, 0, 0.01, 0, 0, 1000],
    "Honey": [6, 100, 100, 95, 0.01, 0, 0, 1000],
    "Jelly": [6.9, 40, 40, 0, 0.01, 7, .7, 1000],
    "Air": [0, 120, 183, 14, 0.005, 0, 0, 0],
};

function selectPreset(preset) {
    document.getElementById("dropdown-button-text").innerHTML = preset;
    document.getElementById("dropdown-button-arrow").innerHTML = "arrow_drop_down";

    /*     for (let i = 0; i < sliders.length; i++) {
            sliders[i].changeStartingValue(presets[preset][i]);
        } */
}

window.onclick = function (event) {
    if (!(event.target.matches('#dropdown-button') || event.target.matches("#dropdown-button-text") || event.target.matches("#dropdown-button-arrow"))) {
        var dropdowns = document.getElementsByClassName("dropdown-content");
        var i;
        for (i = 0; i < dropdowns.length; i++) {
            var openDropdown = dropdowns[i];
            if (openDropdown.classList.contains('dropdown-show')) {
                openDropdown.classList.remove('dropdown-show');
                document.getElementById("dropdown-button-arrow").innerHTML = "arrow_drop_down";
            }
        }
    }
}

function start() {
    running = true;
}

function stopp() {
    running = false;
}

function stepOne() {
    fluid.update(delta);
    fluid.draw();
}

function restart() {
    /* const rD = fluid.restDensity;
    const stiff = fluid.stiffness;
    const stiffN = fluid.nearStiffness;
    const visL = fluid.linearViscosity;
    const visQ = fluid.quadraticViscostiy;
    const springS = fluid.springStiffness;
    const yield = fluid.yieldRate;
    const alpha = fluid.plasticity;
    const gravity = fluid.gravity;

    fluid = new Fluid(fluidAmount);
    fluid.restDensity = rD;
    fluid.stiffness = stiff;
    fluid.nearStiffness = stiffN;
    fluid.linearViscosity = visL;
    fluid.quadraticViscostiy = visQ;
    fluid.springStiffness = springS;
    fluid.yieldRate = yield;
    fluid.plasticity = alpha;
    fluid.gravity = gravity;

    fluid.draw(); */
    fluid = new Fluid(resolutionX, resolutionY, density);
    //fluid.initGrid(10, 10);
    fluid.draw();
}

function test() {
    return [1, 3];
}

function changeResolutionX(newVal) {
    resolutionX = parseInt(newVal);
    fluid = new Fluid(resolutionX, resolutionY, density);
}

function changeResolutionY(newVal) {
    resolutionY = parseInt(newVal);
    fluid = new Fluid(resolutionX, resolutionY, density);
}

function changeDensity(newVal) {
    density = parseFloat(newVal);
    fluid = new Fluid(resolutionX, resolutionY, density);
}

function drawVelocityChangeX(newVal) {
    drawVelX = newVal;
}

function drawVelocityChangeY(newVal) {
    drawVelY = newVal;
}

function drawRedChange(newVal) {
    drawRed = newVal;
}

let frameCount = 0;

function saveCanvasAsPng(filename = "frame.png") {
    const canvas = document.getElementById("canvas");

    canvas.toBlob(function (blob) {
        const a = document.createElement("a");
        document.body.appendChild(a);
        a.style = "display: none";

        const url = window.URL.createObjectURL(blob);
        a.href = url;
        a.download = filename;
        a.click();

        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);
    }, "image/png");
}

function step(now) {
    if (!previous) { previous = now; };

    delta = (now - previous) * 0.001; // seconds
    previous = now;

    if (newFps > 0.5) {
        const average = array => array.reduce((a, b) => a + b) / array.length;

        fps.innerHTML = "" + Math.floor(average(accumulatedFps)) + "fps";
        newFps = 0;
        accumulatedFps = [];
    }
    accumulatedFps.push(1 / delta);
    newFps += delta;

    const minDelta = 0.0069;
    const maxDelta = 0.03;
    delta = Math.max(minDelta, delta);
    delta = Math.min(maxDelta, delta);

    /*     if (delta == minDelta || delta == maxDelta) {
            console.log("delta was adjusted")
        } */

    //running = false;
    delta = FIX_DELTA;

    if (running) {
        fluid.addInflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
        fluid.update(delta);

        fluid.draw(delta);

/*         if (frameCount % 100 === 0) {
            saveCanvasAsPng(`frames/gs${String(frameCount).padStart(6, '0')}.png`);
        } */

        frameCount++;
    }

    requestAnimationFrame(step);
}

requestAnimationFrame(step);