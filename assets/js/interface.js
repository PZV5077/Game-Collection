// Generate dot grid elements
const dotGrid = document.querySelector('.dot-grid');
const dotSize = 20; // Size of each dot in pixels
const dotSpacing = 20; // Spacing between dots in pixels
const dotColor = '#fff'; // Color of the dots
const screenWidth = window.innerWidth;
const screenHeight = window.innerHeight;
const numDotsX = Math.ceil(screenWidth / dotSpacing);
const numDotsY = Math.ceil(screenHeight / dotSpacing);
for (let i = 0; i < numDotsX; i++) {
    for (let j = 0; j < numDotsY; j++) {
        const dot = document.createElement('div');
        dot.classList.add('dot');
        dot.style.top = `${j * dotSpacing}px`;
        dot.style.left = `${i * dotSpacing}px`;
        dot.style.backgroundColor = dotColor;
        dotGrid.appendChild(dot);
    }
}
// Adjust the position of the game button labels
const gameButtons = document.querySelectorAll('.game-button');
gameButtons.forEach(button => {
    const index = button.getAttribute('data-index');
    const label = button.querySelector('.game-button::before');
    label.style.top = `-${label.offsetHeight / 2}px`;
    label.style.left = `-${label.offsetWidth / 2}px`;
});