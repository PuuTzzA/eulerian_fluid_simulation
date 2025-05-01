// import { Fluid } from "./logic.js";
window.addEventListener("keydown", (event) => {

});

const IDLE = 0;
const DRAWING_INFLOW = 1;
const EDITING_INFLOW = 2;

let mode = IDLE;
let x1, y1;
let activeInflow = -1;
const inflowOptionsContainer = document.getElementById("inflow-options");

function getXYfromEvent(e) {
    const x = Math.floor(e.offsetX / fluid.gridPixelSize);
    const offsetTop = window.innerHeight - fluid.gridPixelSize * fluid.h;
    const y = Math.floor((e.offsetY - offsetTop) / fluid.gridPixelSize);
    return [x, y];
}

function openInflowOptions(newInflow) {
    const offsetTop = window.innerHeight - fluid.gridPixelSize * fluid.h;
    const fac = Math.min(fluid.w, fluid.h);
    inflowOptionsContainer.style.top = (activeInflow.y1 * fac * fluid.gridPixelSize + offsetTop) + "px";
    inflowOptionsContainer.style.left = (activeInflow.x1 * fac * fluid.gridPixelSize) + "px";
    inflowOptionsContainer.style.display = "block";

    if (newInflow) {
        inflowOptionsContainer.children[0].resetSlider();
        inflowOptionsContainer.children[1].resetSlider();
        inflowOptionsContainer.children[2].resetSlider();
        inflowOptionsContainer.children[3].resetSlider();
    } else {
        inflowOptionsContainer.children[0].setValue(activeInflow.d);
        inflowOptionsContainer.children[1].setValue(activeInflow.t);
        inflowOptionsContainer.children[2].setValue(activeInflow.u);
        inflowOptionsContainer.children[3].setValue(activeInflow.v);
    }
}

function closeInflowOptions() {
    inflowOptionsContainer.style.display = "none";
}

window.addEventListener("pointerdown", e => {
    if (e.target.id != "canvas") return;
    [x1, y1] = getXYfromEvent(e);
    const fac = Math.min(fluid.w, fluid.h);
    x1 /= fac;
    y1 /= fac;

    if (activeInflow != -1) {
        activeInflow.active = false;
        closeInflowOptions();
    }

    activeInflow = fluid.getInflowAtPoint(x1, y1);
    mode = activeInflow == -1 ? DRAWING_INFLOW : EDITING_INFLOW;

    if (mode == DRAWING_INFLOW) {
        activeInflow = new Inflow(x1, y1, x1, y1, 1, fluid.tAmb, 0, 0, true);
        fluid.inflows.push(activeInflow);
    } else if (mode == EDITING_INFLOW) {
        openInflowOptions(false);
    }
});

window.addEventListener("pointerup", (e) => {
    if (e.target.id != "canvas") return;
    const fac = Math.min(fluid.w, fluid.h);

    if (mode == DRAWING_INFLOW) {
        if (activeInflow.x1 == activeInflow.x2 || activeInflow.y1 == activeInflow.y2) {
            const index = fluid.inflows.indexOf(activeInflow);
            fluid.inflows.splice(index, 1);
            console.log("deleted one inflow");
        } else {
            openInflowOptions(true);
        }
        mode = EDITING_INFLOW;
    } else if (mode == EDITING_INFLOW) {

    }
});


document.addEventListener("mousemove", (e) => {
    if (e.target.id != "canvas") { return };
    const fac = Math.min(fluid.w, fluid.h);

    if (mode == DRAWING_INFLOW) {
        let [x2, y2] = getXYfromEvent(e);

        let _x1 = Math.min(x1, x2 / fac);
        let _x2 = Math.max(x1, x2 / fac);
        let _y1 = Math.min(y1, y2 / fac);
        let _y2 = Math.max(y1, y2 / fac);

        activeInflow.x1 = _x1;
        activeInflow.y1 = _y1;
        activeInflow.x2 = _x2;
        activeInflow.y2 = _y2;
    }
})

function drawInflowD(newVal) {
    activeInflow.d = parseFloat(newVal);
}

function drawInflowT(newVal) {
    activeInflow.t = parseFloat(newVal);
}

function drawInflowU(newVal) {
    activeInflow.u = parseFloat(newVal);
}

function drawInflowV(newVal) {
    activeInflow.v = parseFloat(newVal);
}

window.addEventListener("pointerup", e => {
})

window.addEventListener("resize", e => {
    //fluid.initGrid(fluid.resolutionX, fluid.resolutionY);
})

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
    newFluid();
    //fluid.initGrid(10, 10);
    fluid.draw();
}

function test() {
    return [1, 3];
}

function changeResolutionX(newVal) {
    resolutionX = parseInt(newVal);
    newFluid();
}

function changeResolutionY(newVal) {
    resolutionY = parseInt(newVal);
    newFluid();
}

function changeDensity(newVal) {
    density = parseFloat(newVal);
    newFluid();
}

function changeDensityAir(newVal) {
    densityAir = newVal;
    fluid.densityAir = densityAir;
}

function changeDensitySoot(newVal) {
    densitySoot = newVal;
    fluid.densitySoot = densitySoot;
}

function changeDiffusion(newVal) {
    diffusion = newVal;
    fluid.diffusion = newVal;
}

function changeGravityX(newVal) {
    gravityX = newVal;
    fluid.gravityX = newVal;
}

function changeGravityY(newVal) {
    gravityY = newVal;
    fluid.gravityY = newVal;
}

const fps = document.getElementById("fps");
let newFps = 0;
let accumulatedFps = [];
let running = true;
let delta = 0.0025;
let previous;

let resolutionX = 1;
let resolutionY = 1;
let density = 0;
let densityAir = 0.1;
let densitySoot = 0.1;
let diffusion = 0.01;

let bodies = [];
bodies.push(new SolidBox(0.5, 0.6, 0.5, 0.1, Math.PI * 0.25, 0, 0, 5));
bodies.push(new SolidSphere(0.2, 0.2, 0.2, 0, 0, 0, 0));

let inflows = [];
inflows.push(new Inflow(.38, .1, .62, .2, 1.5, 294, 0.0, 0.0, false));

// bodies.push(new SolidBox(0.5, 0.6, 0.7, 0.1, Math.PI * 0.25, 0.0, 0.0, 0.0));
function newFluid() {
    fluid = new Fluid(resolutionX, resolutionY, densityAir, densitySoot, diffusion, bodies, inflows);
}
let fluid;
newFluid();


const DEFAULT_DT = 0.0025;
const DEFAULT_CFL = 0.5;
const modeInidator = document.getElementById("mode");

function changeDtToCfl() {
    setDt = !setDt;

    if (setDt) {
        document.getElementById("detailed-options").children[8].changeTitle("Timestep");
        document.getElementById("detailed-options").children[8].changeBounds(0.0001, 0.2);
        document.getElementById("detailed-options").children[8].changeStartingValue(DEFAULT_DT);
    } else {
        document.getElementById("detailed-options").children[8].changeTitle("CFL number");
        document.getElementById("detailed-options").children[8].changeBounds(0.001, 20);
        document.getElementById("detailed-options").children[8].changeStartingValue(DEFAULT_CFL);
    }
}

let dtOrCFL = 0.0025;
let specifiedDt = 0.0025;

function changeDTofCFL(newVal) {
    dtOrCFL = newVal;
}


/* fluid.setResolutionX(5);
fluid.setResolutionY(5);
fluid.update(FIX_DELTA); */
fluid.draw();

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
    if (running) {
        fluid.update(specifiedDt);

        bodies.forEach(body => {
            body.update(specifiedDt);
        })

        if (setDt) {
            const cflNumber = (dtOrCFL * fluid.maxVel()) / fluid.cellSize;
            modeInidator.innerHTML = "CFL number: " + cflNumber;
            specifiedDt = dtOrCFL;
        } else {
            const dt = dtOrCFL * fluid.cellSize / fluid.maxVel();
            modeInidator.innerHTML = "Timestep: " + dt;
            specifiedDt = dt;
        }

        /*         if (frameCount % 100 === 0) {
                    saveCanvasAsPng(`frames/gs${String(frameCount).padStart(6, '0')}.png`);
                } */

        frameCount++;
    }

    fluid.draw(delta);


    requestAnimationFrame(step);
}

requestAnimationFrame(step);