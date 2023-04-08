const positionInput = document.getElementById('position');
const addButton = document.getElementById('add');
const removeButton = document.getElementById('remove');
function textFieldChanged(row, textField) {
    var newValue = textField.value;
    google.colab.kernel.invokeFunction("text_callback", [row, newValue], {});
}
function toggleCellContent(cell, row, col) {
  if (row === col) {
    const textCell = document.getElementById(`cell_${row}`);
    const diagMap = { 'H': ['E', 'yellow', 5], 'E': ['C', 'lime', 3], 'C': ['?', 'lightgray', 0], '?': ['H', 'red', 19],  };
    [cell.textContent, cell.style.backgroundColor, textCell.value] = diagMap[cell.textContent] || ['', ''];
    if (cell.textContent === "?" || cell.textContent === "H" || cell.textContent === "C"){
      const gridContainer = document.querySelector('.grid-container');
      const currentSize = (Math.sqrt(1 + 4 * gridContainer.childElementCount) - 1) / 2 - 1;
      for (let k = 0; k < currentSize; k++){
        if (k !== row){
          const a = document.getElementById(`cell_${k}_${row}`);
          const b = document.getElementById(`cell_${row}_${k}`);
          const c = document.getElementById(`cell_${k}_${k}`);
          if (cell.textContent === "?"){
            a.textContent = b.textContent = '?';
            a.style.backgroundColor = b.style.backgroundColor = 'lightgray';
            a.style.opacity = b.style.opacity = 0.1;
          } else if (c.textContent !== "?"){
            a.textContent = b.textContent = '0';
            a.style.backgroundColor = b.style.backgroundColor = 'white';
            if (cell.textContent === "C"){
              a.style.opacity = b.style.opacity = 0.1;
            } else {
              a.style.opacity = b.style.opacity = 1.0;
            }
          }
        }
      }
    }
  } else {
    const diagRow = document.getElementById(`cell_${row}_${row}`);
    const diagCol = document.getElementById(`cell_${col}_${col}`);
    if (diagRow.textContent !== "?" && diagCol.textContent !== "?" && diagRow.textContent !== "C" && diagCol.textContent !== "C"){
      const symcell = document.getElementById(`cell_${col}_${row}`);
      const offdiagMap = { '0': ['1', 'lightblue'], '1': ['?', 'lightgray'], '?': ['0', 'white'] };
      [cell.textContent, cell.style.backgroundColor] = [symcell.textContent, symcell.style.backgroundColor] = offdiagMap[cell.textContent] || ['', ''];
    }
  }
  google.colab.kernel.invokeFunction("toggle_callback", [row, col], {});
}
function updateGridClickEvents() {
  const cells = document.querySelectorAll('.grid-item');
  cells.forEach(cell => {
    cell.addEventListener('click', () => {
      const [row, col] = cell.id.split('_').slice(1).map(Number);
      toggleCellContent(cell, row, col);
    });
  });
}
function createGridItem(row, col, bgColor, content, opacity) {
  return `<div class="grid-item" id="cell_${row}_${col}" style="background-color:${bgColor};opacity:${opacity}">${content}</div>`;
}
function createTextInput(row, value) {
  return `<div><input class="text" type="number" id="cell_${row}" min="0" value="${value}" onchange="textFieldChanged(${row}, this)"></input></div>`;
}
function updateGrid(add) {
  const gridContainer = document.querySelector('.grid-container');
  const currentSize = (Math.sqrt(1 + 4 * gridContainer.childElementCount) - 1) / 2 - 1;
  const newSize = add ? currentSize + 1 : currentSize - 1;
  let position = parseInt(positionInput.value);
  if (position === -1) { position = currentSize; }

  let newGrid = '<div class="pos"></div>';
  for (let row = 0; row < newSize; row++) {
    newGrid += `<div class="pos">${row}</div>`;
  }
  newGrid += '<div class="pos"></div>';
  for (let row = 0; row < newSize; row++) {
    const oldRow = add ? (row < position ? row : row - 1) : (row < position ? row : row + 1);
    newGrid += `<div class="pos">${row}</div>`;
    for (let col = 0; col < newSize; col++) {
      const oldCol = add ? (col < position ? col : col - 1) : (col < position ? col : col + 1);
      let bgColor, content, opacity;
      if (add && (row === position || col === position)) {
        if (row === col) {
          bgColor = 'red';
          content = 'H';
          opacity = 1.0;
        } else {
          let cell = col === position ? document.getElementById(`cell_${oldRow}_${oldRow}`) : document.getElementById(`cell_${oldCol}_${oldCol}`);
          if (cell.textContent === "?") {
            bgColor = 'lightgray';
            content = '?';
            opacity = 0.1;
          } else {
            bgColor = 'white';
            content = '0';
            opacity = cell.textContent === "C" ? 0.1 : 1.0;
          }
        }
        newGrid += createGridItem(row, col, bgColor, content, opacity);
      } else {
        const cell = document.getElementById(`cell_${oldRow}_${oldCol}`);
        newGrid += createGridItem(row, col, cell.style.backgroundColor, cell.textContent, cell.style.opacity);
      }
    }
    if (add && row === position) {
      newGrid += createTextInput(row, 19);
    } else {
      const cell = document.getElementById(`cell_${oldRow}`);
      newGrid += createTextInput(row, cell.value);
    }
  }
  gridContainer.style.gridTemplateColumns = `repeat(${newSize + 2}, 30px)`;
  gridContainer.innerHTML = newGrid;
  updateGridClickEvents();
  positionInput.max = newSize;
  google.colab.kernel.invokeFunction("update_callback", [position, add], {});
}
addButton.addEventListener('click', function()    {updateGrid(true);});
removeButton.addEventListener('click', function() {updateGrid(false);});
updateGridClickEvents();