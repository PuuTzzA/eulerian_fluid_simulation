const canvas = document.getElementById("canvas");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
const ctx = canvas.getContext("2d");

const FLUID = 1;
const SOLID = 2;
const EMPTY = 3;

window.addEventListener("resize", e => {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
})

class GraphicsSettings {
    constructor() {
        this.lineWidth = 3;
        this.lineColour = "white";
        this.particleSize = 3;
        this.maxSpeed = 10;
        this.maxVelocitySize = 1;
        this.arrowWidth = 2;
        this.arrowHeadSize = 9;
        this.arrowColor = "red";
    }

    getLineWidth() {
        return this.lineWidth;
    }

    getLineColour() {
        return this.lineColour;
    }

    getHatchingSettings(label, gridSize) {
        const solidSettings = {
            colour: "rgb(112, 101, 84)",
            lineWidth: 3,
            cross: true,
            direction: "right",
            delta: gridSize / 6
        };
        const fluidSettings = {
            colour: "rgb(28, 164, 255)",
            lineWidth: 1,
            cross: false,
            direction: "left",
            delta: gridSize / 4
        };
        const emptySettings = {
            colour: "white",
            lineWidth: 0.5,
            cross: false,
            direction: "right",
            delta: gridSize / 2
        };

        switch (label) {
            case FLUID:
                return fluidSettings;
            case SOLID:
                return solidSettings;
            case EMPTY:
                return emptySettings;
            default:
                return "rgb(255, 0, 255)";
        }
    }

    getHatchingDelta() {

    }

    getParticleRadius() {
        return this.particleSize;
    }

    getCellColour(red, green, blue) {
        return `rgb(${red}, ${green}, ${blue})`;
    }

    getArrowWidth() {
        return this.arrowWidth;
    }

    getArrowHeadSize() {
        return this.arrowHeadSize;
    }

    getArrowColor(velocity) {
        return this.arrowColor;
    }

    getArrowLength(velocity) {
        // returns a fraction of the grid spacing
        return (velocity / this.maxSpeed) * this.maxVelocitySize;
    }
}

class Fluid {
    constructor(resolutionX, resolutionY) {
        this.gs = new GraphicsSettings();

        this.resolutionX = resolutionX;
        this.resolutionY = resolutionY;

        this.mousePos = [0, 0];
        this.mousePosPrev = [0, 0];
        this.mouseRadius = 100;
        this.mouseRadiusSquared = this.mouseRadius * this.mouseRadius;
        this.mousePressed = false;
        this.mouseForce = 200;

        this.gravity = 2000;

        canvas.width = window.innerWidth;
        canvas.height = window.innerHeight;

        this.density = 1000;
        this.densityReciprocal = 1 / this.density;

        this.pressures = new Float32Array(1);
        this.levelSet = new Float32Array(1);
        this.labels = new Int8Array(1);
        this.negativeDivergence = new Float32Array(1);
        this.Adiag = new Float32Array(1);
        this.Ax = new Float32Array(1);
        this.Ay = new Float32Array(1);
        this.preconditioner = new Float64Array(1);
        this.velocitiesX = new Float32Array(1);
        this.velocitiesY = new Float32Array(1);
        this.colorRed = new Float32Array(1);

        this.gridSize = 10;
    }

    setResolutionX(newX) {
        this.initGrid(parseInt(newX), this.resolutionY);
    }

    setResolutionY(newY) {
        this.initGrid(this.resolutionX, parseInt(newY));
    }

    initGrid(resX, resY) {
        let shouldX = window.innerWidth / resX;
        let shouldY = window.innerHeight / resY;

        this.gridSize = Math.min(shouldX, shouldY);
        this.resolutionX = resX;
        this.resolutionY = resY;

        this.pressures = new Float32Array(resX * resY);
        this.levelSet = new Float32Array(resX * resY);
        this.labels = new Int8Array(resX * resY);
        this.negativeDivergence = new Float64Array(resX * resY);
        this.Adiag = new Float64Array(resX * resY);
        this.Ax = new Float64Array(resX * resY);
        this.Ay = new Float64Array(resX * resY);
        this.preconditioner = new Float64Array(resX * resY);
        this.colorRed = new Float32Array(resX * resY);
        this.velocitiesX = new Float32Array((resX + 1) * resY);
        this.velocitiesY = new Float32Array(resX * (resY + 1));

        this.velocitiesX[(resX + 1) * 4 + 3] = 2;
        this.velocitiesX[(resX + 1) * 5 + 3] = 2;
        this.velocitiesX[(resX + 1) * 6 + 3] = 2;

        for (let i = 0, n = this.labels.length; i < n; i++) {
            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            if (x == 0 || x == this.resolutionX - 1 || y == this.resolutionY - 1) {
                this.labels[i] = SOLID;
            }

            this.levelSet[i] = Fluid.aaSdf(2, resX - 2, 1, resY - 3, x, y);
        }
    }

    static aaSdf(x0, x1, y0, y1, x, y) {
        if (x0 < x && x < x1 && y0 < y && y < y1) {
            // inside the box
            return Math.max(x0 - x, x - x1, y0 - y, y - y1);
        } else {
            // outside
            let p, q;
            if (x < x0) { p = x0; }
            else if (x > x1) { p = x1; }
            else { p = x; };
            if (y < y0) { q = y0; }
            else if (y > y1) { q = y1; }
            else { q = y; };
            return Math.sqrt((x - p) * (x - p) + (y - q) * (y - q));
        }
    }

    getVelocityAtPoint(x, y) {
        const u = Fluid.bilinearInterpolation(x + 0.5, y, this.velocitiesX, this.resolutionX + 1, this.resolutionY);
        const v = Fluid.bilinearInterpolation(x, y + 0.5, this.velocitiesY, this.resolutionX, this.resolutionY + 1);
        return [u, v];
    }

    static bilinearInterpolation(x, y, values, resolutionX, resolutionY) {
        // Clamp to avoid out-of-bounds access
        x = Math.max(0, Math.min(x, resolutionX - 1));
        y = Math.max(0, Math.min(y, resolutionY - 1));

        const x0 = Math.floor(x);
        const y0 = Math.floor(y);
        const x1 = Math.min(x0 + 1, resolutionX - 1);
        const y1 = Math.min(y0 + 1, resolutionY - 1);

        const alpha = x - x0;
        const beta = y - y0;

        const topLeft = values[y0 * resolutionX + x0];
        const topRight = values[y0 * resolutionX + x1];
        const bottomLeft = values[y1 * resolutionX + x0];
        const bottomRight = values[y1 * resolutionX + x1];

        const top = topLeft * (1 - alpha) + topRight * alpha;
        const bottom = bottomLeft * (1 - alpha) + bottomRight * alpha;

        return top * (1 - beta) + bottom * beta;
    }

    static bicubicInterpolation(x, y, values, resolutionX, resolutionY) {
        x = Math.max(0, Math.min(x, resolutionX - 1));
        y = Math.max(0, Math.min(y, resolutionY - 1));

        const x1 = Math.floor(x);
        const y1 = Math.floor(y);

        const x0 = Math.max(x1 - 1, 0);
        const x2 = Math.min(x1 + 1, resolutionX - 1);
        const x3 = Math.min(x1 + 2, resolutionX - 1);

        const y0 = Math.max(y1 - 1, 0);
        const y2 = Math.min(y1 + 1, resolutionY - 1);
        const y3 = Math.min(y1 + 2, resolutionY - 1);

        const alpha = x - x1;
        const beta = y - y1;

        // weighing coefficients
        const w0 = (-1 / 3) * alpha + 0.5 * alpha * alpha - (1 / 6) * alpha * alpha * alpha;
        const w1 = 1 - alpha * alpha + 0.5 * (alpha * alpha * alpha - alpha);
        const w2 = alpha + 0.5 * (alpha * alpha - alpha * alpha * alpha);
        const w3 = (1 / 6) * (alpha * alpha * alpha - alpha);

        const ww0 = (-1 / 3) * beta + 0.5 * beta * beta - (1 / 6) * beta * beta * beta;
        const ww1 = 1 - beta * beta + 0.5 * (beta * beta * beta - beta);
        const ww2 = beta + 0.5 * (beta * beta - beta * beta * beta);
        const ww3 = (1 / 6) * (beta * beta * beta - beta);

        // top top row (1st)
        const topTopLeftLeft = values[y0 * resolutionX + x0];
        const topTopLeft = values[y0 * resolutionX + x1];
        const topTopRight = values[y0 * resolutionX + x2];
        const topTopRightRight = values[y0 * resolutionX + x3];

        const topTop = w0 * topTopLeftLeft + w1 * topTopLeft + w2 * topTopRight + w3 * topTopRightRight;

        // top row (2nd)
        const topLeftLeft = values[y1 * resolutionX + x0];
        const topLeft = values[y1 * resolutionX + x1];
        const topRight = values[y1 * resolutionX + x2];
        const topRightRight = values[y1 * resolutionX + x3];

        const top = w0 * topLeftLeft + w1 * topLeft + w2 * topRight + w3 * topRightRight;

        // bottom row (3rd)
        const bottomLeftLeft = values[y2 * resolutionX + x0];
        const bottomLeft = values[y2 * resolutionX + x1];
        const bottomRight = values[y2 * resolutionX + x2];
        const bottomRightRight = values[y2 * resolutionX + x3];

        const bottom = w0 * bottomLeftLeft + w1 * bottomLeft + w2 * bottomRight + w3 * bottomRightRight;

        // bottom bottom row (4th)
        const bottomBottomLeftLeft = values[y3 * resolutionX + x0];
        const bottomBottomLeft = values[y3 * resolutionX + x1];
        const bottomBottomRight = values[y3 * resolutionX + x2];
        const bottomBottomRightRight = values[y3 * resolutionX + x3];

        const bottomBottom = w0 * bottomBottomLeftLeft + w1 * bottomBottomLeft + w2 * bottomBottomRight + w3 * bottomBottomRightRight;

        // interpolation of the rows
        return ww0 * topTop + ww1 * top + ww2 * bottom + ww3 * bottomBottom;
    }

    advect(dt, oldVals, resolutionX, resolutionY) {
        let newVals = new Float32Array(resolutionX * resolutionY);

        for (let i = 0, n = newVals.length; i < n; i++) {
            // if (i != 0) { continue };
            // Runge Kutta 4
            const x = i % resolutionX;
            const y = Math.floor(i / resolutionX);

            const k1 = this.getVelocityAtPoint(x, y);
            const k2 = this.getVelocityAtPoint(x - k1[0] * dt / 2, y - k1[1] * dt / 2);
            const k3 = this.getVelocityAtPoint(x - k2[0] * dt / 2, y - k2[1] * dt / 2);
            const k4 = this.getVelocityAtPoint(x - k3[0] * dt, y - k3[1] * dt);

            let newX = x - (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
            let newY = y - (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);

            newX = Math.max(newX, 0);
            newX = Math.min(newX, resolutionX - 1);
            newY = Math.max(newY, 0);
            newY = Math.min(newY, resolutionY - 1);

            newVals[i] = Fluid.bicubicInterpolation(newX, newY, oldVals, resolutionX, resolutionY);
        }

        return newVals;
    }

    updateLabels() {
        for (let i = 0, n = this.labels.length; i < n; i++) {
            if (this.labels[i] == SOLID) { continue; }

            if (this.levelSet[i] <= 0) {
                this.labels[i] = FLUID;
            } else if (this.levelSet[i] > 0) {
                this.labels[i] = EMPTY;
            }
        }
    }

    fastSweep2D(distGrid, width, height) {
        const NSweeps = 4;
        const h = 1.0;
        const f = 1.0;
        const eps = 1e-6;

        const dirX = [
            [0, width - 1, 1],
            [width - 1, 0, -1],
            [width - 1, 0, -1],
            [0, width - 1, 1]
        ];
        const dirY = [
            [0, height - 1, 1],
            [0, height - 1, 1],
            [height - 1, 0, -1],
            [height - 1, 0, -1]
        ];

        for (let s = 0; s < NSweeps; s++) {
            const [x0, x1, xStep] = dirX[s];
            const [y0, y1, yStep] = dirY[s];

            for (let iy = y0; yStep * iy <= yStep * y1; iy += yStep) {
                for (let ix = x0; xStep * ix <= xStep * x1; ix += xStep) {
                    const gridPos = iy * width + ix;

                    // Find min neighbor in x-direction
                    let ax;
                    if (ix === 0) {
                        ax = Math.min(distGrid[gridPos], distGrid[iy * width + (ix + 1)]);
                    } else if (ix === width - 1) {
                        ax = Math.min(distGrid[iy * width + (ix - 1)], distGrid[gridPos]);
                    } else {
                        ax = Math.min(distGrid[iy * width + (ix - 1)], distGrid[iy * width + (ix + 1)]);
                    }

                    // Find min neighbor in y-direction
                    let ay;
                    if (iy === 0) {
                        ay = Math.min(distGrid[gridPos], distGrid[(iy + 1) * width + ix]);
                    } else if (iy === height - 1) {
                        ay = Math.min(distGrid[(iy - 1) * width + ix], distGrid[gridPos]);
                    } else {
                        ay = Math.min(distGrid[(iy - 1) * width + ix], distGrid[(iy + 1) * width + ix]);
                    }

                    // Solve local Eikonal update
                    let d_new;
                    if (Math.abs(ax - ay) < f * h) {
                        d_new = 0.5 * (ax + ay + Math.sqrt(2 * f * f * h * h - (ax - ay) * (ax - ay)));
                    } else {
                        d_new = Math.min(ax, ay) + f * h;
                    }

                    // Only accept the update if it improves the value
                    distGrid[gridPos] = Math.min(distGrid[gridPos], d_new);
                }
            }
        }
    }

    applyBodyForces(dt) {
        for (let i = 0, n = this.velocitiesY.length; i < n; i++) {
            this.velocitiesY[i] += this.gravity * dt;
        }
    }

    pressureGradientUpdate(dt) {
        const scale = dt / (this.density * 1) // 1 == dx

        for (let i = 0; i < this.pressures.length; i++) {
            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            // update u
            if (x > 0 &&
                (this.labels[y * this.resolutionX + x - 1] == FLUID || this.labels[y * this.resolutionX + x] == FLUID)) {

                if (this.labels[y * this.resolutionX + x - 1] == SOLID || this.labels[y * this.resolutionX + x] == SOLID) {
                    this.velocitiesX[y * (this.resolutionX + 1) + x] = 0; // u(i, j) = u_solid(i, j)
                } else {
                    const p_left = this.pressures[y * this.resolutionX + x - 1];
                    const p_right = this.pressures[y * this.resolutionX + x];
                    this.velocitiesX[y * (this.resolutionX + 1) + x] -= scale * (p_right - p_left);
                }
            } else {
                // this.velocitiesX[y * this.resolutionX + x] = UNKNOWN;
                this.velocitiesX[y * (this.resolutionX + 1) + x] = 0;
            }

            // update v
            if (y > 0 &&
                (this.labels[(y - 1) * this.resolutionX + x] == FLUID || this.labels[y * this.resolutionX + x] == FLUID)) {

                if (this.labels[(y - 1) * this.resolutionX + x] == SOLID || this.labels[y * this.resolutionX + x] == SOLID) {
                    this.velocitiesY[y * this.resolutionX + x] = 0; // u(i, j) = u_solid(i, j)
                } else {
                    const p_top = this.pressures[(y - 1) * this.resolutionX + x];
                    const p_bottom = this.pressures[y * this.resolutionX + x];
                    this.velocitiesY[y * this.resolutionX + x] -= scale * (p_bottom - p_top);
                }
            } else {
                // this.velocitiesY[y * this.resolutionX + x] = UNKNOWN;
                this.velocitiesY[i] = 0;
            }
        }
    }

    calculateNegativeDivergence() {
        const scale = 1 / 1; // 1/dx

        for (let i = 0, n = this.negativeDivergence.length; i < n; i++) {
            if (this.labels[i] != FLUID) { continue; }

            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            const du = this.velocitiesX[y * (this.resolutionX + 1) + x + 1] - this.velocitiesX[y * (this.resolutionX + 1) + x];
            const dv = this.velocitiesY[(y + 1) * this.resolutionX + x] - this.velocitiesY[y * this.resolutionX + x];

            this.negativeDivergence[i] = -scale * (du + dv);

            if (x > 0 && this.labels[y * this.resolutionX + x - 1] == SOLID) {
                this.negativeDivergence[i] -= scale * this.velocitiesX[y * (this.resolutionX + 1) + x] - 0; // u(x, y) - u_solid(x, y)
            }
            if (x < this.resolutionX - 1 && this.labels[y * this.resolutionX + x + 1] == SOLID) {
                this.negativeDivergence[i] += scale * this.velocitiesX[y * (this.resolutionX + 1) + x + 1] - 0; // u(x, y) - u_solid(x, y)
            }

            if (y > 0 && this.labels[(y - 1) * this.resolutionX + x] == SOLID) {
                this.negativeDivergence[i] -= scale * this.velocitiesY[y * this.resolutionX + x] - 0; // v(x, y) - v_solid(x, y)
            }
            if (y < this.resolutionY - 1 && this.labels[(y + 1) * this.resolutionX + x] == SOLID) {
                this.negativeDivergence[i] += scale * this.velocitiesY[(y + 1) * this.resolutionX + x] - 0; // v(x, y) - v_solid(x, y)
            }
        }
    }

    setupA(dt) {
        const scale = dt / (this.density * 1 * 1); // dt / (density * dx^2)
        this.Adiag.fill(0);
        this.Ax.fill(0);
        this.Ay.fill(0);

        for (let i = 0, n = this.pressures.length; i < n; i++) {
            if (this.labels[i] == FLUID) {
                const x = i % this.resolutionX;
                const y = Math.floor(i / this.resolutionX);

                // negative x neighbour
                if (x > 0 && this.labels[y * this.resolutionX + x - 1] == FLUID) {
                    this.Adiag[i] += scale;
                }
                else if (x > 0 && this.labels[y * this.resolutionX + x - 1] == EMPTY) {
                    this.Adiag[i] += scale;
                }
                // positive x neighbour
                if (x < this.resolutionX - 1 && this.labels[y * this.resolutionX + x + 1] == FLUID) {
                    this.Adiag[i] += scale;
                    this.Ax[i] = -scale;
                }
                else if (x < this.resolutionX - 1 && this.labels[y * this.resolutionX + x + 1] == EMPTY) {
                    this.Adiag[i] += scale;
                }
                // negative y neighbour
                if (y > 0 && this.labels[(y - 1) * this.resolutionX + x] == FLUID) {
                    this.Adiag[i] += scale;
                }
                else if (y > 0 && this.labels[(y - 1) * this.resolutionX + x] == EMPTY) {
                    this.Adiag[i] += scale;
                }
                // positive y neighbour
                if (y < this.resolutionY - 1 && this.labels[(y + 1) * this.resolutionX + x] == FLUID) {
                    this.Adiag[i] += scale;
                    this.Ay[i] = -scale;
                }
                else if (y < this.resolutionY - 1 && this.labels[(y + 1) * this.resolutionX + x] == EMPTY) {
                    this.Adiag[i] += scale;
                }
            }
        }
    }

    /**
    * The infinity norm is the maximum absolute value of a vector.
    */
    static infinityNorm(s) {
        let max = 0;

        for (let i = 0, n = s.length; i < n; i++) {
            max = Math.max(max, Math.abs(s[i]));
        }

        return max;
    }

    static dotProduct(x, y) {
        let result = 0;

        for (let i = 0, n = x.length; i < n; i++) {
            result += x[i] * y[i];
        }

        return result;
    }

    static addVectors(x, y) {
        let result = new Float64Array(x.length);

        for (let i = 0, n = x.length; i < n; i++) {
            result[i] = x[i] + y[i];
        }

        return result;
    }

    static subtractVectors(x, y) {
        let result = new Float64Array(x.length);

        for (let i = 0, n = x.length; i < n; i++) {
            result[i] = x[i] - y[i];
        }

        return result;
    }

    static multiplyScalarVector(a, x) {
        let result = new Float64Array(x.length);

        for (let i = 0, n = x.length; i < n; i++) {
            result[i] = a * x[i];
        }

        return result;
    }

    /**
    * multiplies the matrix A with the vector s
    */
    applyA(s) {
        let result = new Float64Array(s.length);

        for (let i = 0, n = s.length; i < n; i++) {
            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            let r = this.Adiag[i] * s[i];

            // neighbour at (x - 1, y)
            if (x > 0) {
                r += this.Ax[y * this.resolutionX + x - 1] * s[y * this.resolutionX + x - 1];
            }
            // neighbour at (x + 1, y)
            if (x < this.resolutionX - 1) {
                r += this.Ax[i] * s[y * this.resolutionX + x + 1];
            }
            // neighbour at (x, y - 1)
            if (y > 0) {
                r += this.Ay[(y - 1) * this.resolutionX + x] * s[(y - 1) * this.resolutionX + x];
            }
            // neighbour at (x, y + 1)
            if (y < this.resolutionY - 1) {
                r += this.Ay[i] * s[(y + 1) * this.resolutionX + x];
            }

            result[i] = r;
        }

        return result;
    }

    constructPreconditioner() {
        const tuningConstant = 0.97;
        const safetyConstant = 0.25;

        for (let i = 0, n = this.preconditioner.length; i < n; i++) {
            if (this.labels[i] != FLUID) { continue; }

            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            let e = this.Adiag[i];

            if (x > 0) {
                const s = this.Ax[y * this.resolutionX + x - 1] * this.preconditioner[y * this.resolutionX + x - 1];
                e -= s * s;
            }
            if (y > 0) {
                const s = this.Ay[(y - 1) * this.resolutionX + x] * this.preconditioner[(y - 1) * this.resolutionX + x];
                e -= s * s;
            }
            if (x > 0 && y > 0) {
                const t = this.Ax[y * this.resolutionX + x - 1] * this.Ay[y * this.resolutionX + x - 1]
                    * this.preconditioner[y * this.resolutionX + x - 1] * this.preconditioner[y * this.resolutionX + x - 1]
                    + this.Ay[(y - 1) * this.resolutionX + x] * this.Ax[(y - 1) * this.resolutionX + x]
                    * this.preconditioner[(y - 1) * this.resolutionX + x] * this.preconditioner[(y - 1) * this.resolutionX + x];

                e -= tuningConstant * t;
            }

            if (e < safetyConstant * this.Adiag[i]) {
                e = this.Adiag[i];
            }
            if (e <= 0) { e = 1e-6; }
            this.preconditioner[i] = 1 / Math.sqrt(e);
        }
    }

    applyPreconditioner(r) {
        // first Lq = r  
        let q = new Float64Array(r.length);

        for (let i = 0, n = this.preconditioner.length; i < n; i++) {
            if (this.labels[i] != FLUID) { continue; }

            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            let t = r[i];

            if (x > 0) {
                const index = y * this.resolutionX + x - 1;
                t -= this.Ax[index] * this.preconditioner[index] * q[index];
            }
            if (y > 0) {
                const index = (y - 1) * this.resolutionX + x;
                t -= this.Ay[index] * this.preconditioner[index] * q[index];
            }
            q[i] = t * this.preconditioner[i];
        }
        // solve L^T z = q
        let z = new Float64Array(q.length);

        for (let i = this.preconditioner.length - 1; i >= 0; i--) {
            if (this.labels[i] != FLUID) { continue; }

            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            let t = q[i];

            if (x < this.resolutionX - 1) {
                t -= this.Ax[i] * this.preconditioner[i] * z[y * this.resolutionX + x + 1];
            }
            if (y < this.resolutionY - 1) {
                t -= this.Ay[i] * this.preconditioner[i] * z[(y + 1) * this.resolutionX + x];
            }
            z[i] = t * this.preconditioner[i];
        }
        return z;
    }

    /**
    * The preconditioned conjugate gradient (PCG) algorithm for solving Ap = b.
    */
    preconditionedConjugateGradient() {
        let p = new Float64Array(this.negativeDivergence.length);

        let residual = new Float64Array(this.negativeDivergence);

        if (Fluid.infinityNorm(residual) == 0) {
            return p;
        }

        let auxiliaryVector = this.applyPreconditioner(residual);
        let searchVector = new Float64Array(auxiliaryVector);


        let sigma = Fluid.dotProduct(residual, auxiliaryVector);

        const maxIter = 200;
        const tolerance = 1e-6;

        let alpha = 1;
        let beta = 1;

        for (let _ = 0; _ < maxIter; _++) {
            auxiliaryVector = this.applyA(searchVector);

            alpha = sigma / Fluid.dotProduct(auxiliaryVector, searchVector);

            p = Fluid.addVectors(p, Fluid.multiplyScalarVector(alpha, searchVector));
            residual = Fluid.subtractVectors(residual, Fluid.multiplyScalarVector(alpha, auxiliaryVector));

            if (Fluid.infinityNorm(residual) < tolerance) {
                return p;
            }

            auxiliaryVector = this.applyPreconditioner(residual);

            const sigmaNew = Fluid.dotProduct(auxiliaryVector, residual);
            beta = sigmaNew / sigma;
            searchVector = Fluid.addVectors(auxiliaryVector, Fluid.multiplyScalarVector(beta, searchVector));
            sigma = sigmaNew;
        }

        console.log("PCG LIMIT EXCEEDED!-----------------------------------");
        return p;
    }

    project(dt) {
        this.updateLabels();

        this.calculateNegativeDivergence();
        console.log("divergence", this.negativeDivergence)

        this.setupA(dt);
        console.log("A diag", this.Adiag);
        console.log("A x", this.Ax);
        console.log("A y", this.Ay);

        this.constructPreconditioner();

        this.pressures = this.preconditionedConjugateGradient();

        this.pressureGradientUpdate(dt);
    }

    update(dt) {


        // add u = u + dt * g
        this.applyBodyForces(dt);

        // set u = project(dt, u)
        this.project(dt);

        this.colorRed = this.advect(dt, this.colorRed, this.resolutionX, this.resolutionY);
        this.levelSet = this.advect(dt, this.levelSet, this.resolutionX, this.resolutionY);
        this.fastSweep2D(this.levelSet, this.resolutionX, this.resolutionY);

        // Set u = advect(u, dt)
        const newVelX = this.advect(dt, this.velocitiesX, this.resolutionX + 1, this.resolutionY);
        const newVelY = this.advect(dt, this.velocitiesY, this.resolutionX, this.resolutionY + 1);
        this.velocitiesX = newVelX;
        this.velocitiesY = newVelY;
    }

    draw() {
        ctx.clearRect(0, 0, canvas.width, canvas.height);
        ctx.fillStyle = "black";
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        ctx.lineWidth = this.gs.getLineWidth();

        const offsetTop = window.innerHeight - this.gridSize * this.resolutionY;


        for (let i = 0, n = this.pressures.length; i < n; i++) {
            const x = i % this.resolutionX;
            const y = Math.floor(i / this.resolutionX);

            ctx.strokeStyle = this.gs.getLineColour();
            ctx.fillStyle = this.gs.getCellColour(this.colorRed[i], 0, 0);

            ctx.strokeRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize);
            ctx.fillRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize);

            Fluid.hatchRect(x * this.gridSize, y * this.gridSize + offsetTop, this.gridSize, this.gridSize, this.gs.getHatchingSettings(this.labels[i], this.gridSize));
        }

        // horizontal velocities
        for (let i = 0, n = this.velocitiesX.length; i < n; i++) {
            if (Math.abs(this.velocitiesX[i]) < 0.01) continue;

            const x = i % (this.resolutionX + 1);
            const y = Math.floor(i / (this.resolutionX + 1));

            let fromX = x * this.gridSize;
            let fromY = (y + 0.5) * this.gridSize + offsetTop;
            let toX = (x + this.gs.getArrowLength(this.velocitiesX[i])) * this.gridSize;
            let toY = fromY;

            let col = i == 0 ? "white" : this.gs.getArrowColor();
            Fluid.drawArrow(ctx, fromX, fromY, toX, toY, this.gs.getArrowWidth(), this.gs.getArrowHeadSize(), col);
        }

        // vertical velocities
        for (let i = 0, n = this.velocitiesY.length; i < n; i++) {
            if (Math.abs(this.velocitiesY[i]) < 0.01) continue;

            const x = i % (this.resolutionX);
            const y = Math.floor(i / (this.resolutionX));

            let fromX = (x + 0.5) * this.gridSize;
            let fromY = y * this.gridSize + offsetTop;
            let toX = fromX;
            let toY = (y + this.gs.getArrowLength(this.velocitiesY[i])) * this.gridSize + offsetTop;

            let col = i == 0 ? "white" : this.gs.getArrowColor();

            Fluid.drawArrow(ctx, fromX, fromY, toX, toY, this.gs.getArrowWidth(), this.gs.getArrowHeadSize(), col);
        }
    }

    static drawArrow(ctx, fromx, fromy, tox, toy, arrowWidth, headlen, color) {
        //variables to be used when creating the arrow
        var angle = Math.atan2(toy - fromy, tox - fromx);

        ctx.save();
        ctx.strokeStyle = color;

        //starting path of the arrow from the start square to the end square
        //and drawing the stroke
        ctx.beginPath();
        ctx.moveTo(fromx, fromy);
        ctx.lineTo(tox, toy);
        ctx.lineWidth = arrowWidth;
        ctx.stroke();

        //starting a new path from the head of the arrow to one of the sides of
        //the point
        ctx.beginPath();
        ctx.moveTo(tox, toy);
        ctx.lineTo(tox - headlen * Math.cos(angle - Math.PI / 7),
            toy - headlen * Math.sin(angle - Math.PI / 7));

        //path from the side point of the arrow, to the other side point
        ctx.lineTo(tox - headlen * Math.cos(angle + Math.PI / 7),
            toy - headlen * Math.sin(angle + Math.PI / 7));

        //path from the side point back to the tip of the arrow, and then
        //again to the opposite side point
        ctx.lineTo(tox, toy);
        ctx.lineTo(tox - headlen * Math.cos(angle - Math.PI / 7),
            toy - headlen * Math.sin(angle - Math.PI / 7));

        //draws the paths created above
        ctx.stroke();
        ctx.restore();
    }

    static hatchRect(x1, y1, dx, dy, settings) {
        const cross = settings.cross;
        const direction = settings.direction;
        const color = settings.colour;
        const delta = settings.delta;
        const lw = settings.lineWidth;

        ctx.save();
        ctx.beginPath();
        ctx.rect(x1, y1, dx, dy);
        ctx.clip();

        const majorAxe = Math.max(dx, dy);
        ctx.strokeStyle = color;
        ctx.lineWidth = lw;

        const drawLines = (dir) => {
            for (let n = -majorAxe; n < majorAxe; n += delta) {
                ctx.beginPath();
                if (dir === 'right') {
                    // bottom-left to top-right (like your original)
                    ctx.moveTo(n + x1, y1);
                    ctx.lineTo(n + x1 + dy, y1 + dy);
                } else if (dir === 'left') {
                    // top-left to bottom-right
                    ctx.moveTo(n + x1, y1 + dy);
                    ctx.lineTo(n + x1 + dy, y1);
                }
                ctx.stroke();
            }
        };

        // Draw primary direction
        drawLines(direction);

        // Optionally add crosshatching
        if (cross) {
            drawLines(direction === 'right' ? 'left' : 'right');
        }

        ctx.restore();
    }
}
