// import { Fluid } from "./logic.js";
document.addEventListener("mousemove", (e) => {
})

window.addEventListener("pointerdown", e => {
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
let fluid = new Fluid(1, 1);
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
    // fluid.update(delta);
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
}

function changeResolutionX(newVal) {
    fluid.setResolutionX(newVal);
}

function changeResolutionY(newVal) {
    fluid.setResolutionY(newVal);
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

    if (running) {
        fluid.update(delta);
        fluid.draw(delta);
    }

    requestAnimationFrame(step);
}

requestAnimationFrame(step);