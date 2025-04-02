class MyInputContainer extends HTMLElement {
    connectedCallback() {
        this.innerHTML =
            `
            <div class="container-input">
                ${this.innerHTML}
            </div> 
            `;
    }
}

customElements.define("my-input-container", MyInputContainer);

let sliders = [];

class MySlider extends HTMLElement {
    connectedCallback() {
        const title = this.getAttribute("title");
        const min = this.getAttribute("min");
        const max = this.getAttribute("max");
        const step = this.getAttribute("step");
        this.startingValue = this.getAttribute("startingValue");

        this.innerHTML =
            `
            <div class="slider">
                <div>
                    <span>${title}</span><span class="slider-value">${this.startingValue}</span>
                </div>
                <div>
                    <span class="material-symbols-outlined reset-slider">restart_alt</span>
                    <input class="slider-input" type="range" min="${min}" max="${max}" step="${step}" value="${this.startingValue}">
                </div>
            </div>
            `;

        console.log()

        sliders.push(this)

        this.callbackFunction = this.getAttribute("callbackFunction");
        this.output = this.children[0].children[0].children[1];
        this.slider = this.children[0].children[1].children[1];
        window[this.callbackFunction](this.startingValue);

        this.addEventListener("input", e => {
            this.output.innerHTML = e.target.value;
            window[this.callbackFunction](e.target.value);
        });

        this.children[0].children[1].children[0].addEventListener("click", e => {
            this.output.innerHTML = this.startingValue;
            this.slider.value = this.startingValue;
            window[this.callbackFunction](this.startingValue);
        })
    }

    changeStartingValue(newStartingValue) {
        this.startingValue = newStartingValue;
        this.output.innerHTML = newStartingValue;
        this.slider.value = newStartingValue;
        window[this.callbackFunction](this.startingValue);
    }
}

customElements.define("my-slider", MySlider);

class MyButton extends HTMLElement {
    connectedCallback() {
        const title = this.getAttribute("title");

        this.innerHTML =
            `
                <button class="button material-symbols-outlined" id="button_1">${title}</button>
            `;

        this.callbackFunction = this.getAttribute("callbackFunction");

        this.children[0].addEventListener("click", e => {
            window[this.callbackFunction]("button");
        })
    }
}

customElements.define("my-button", MyButton);