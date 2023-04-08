function textFieldChanged(row, textField) {
    var newValue = textField.value;
    google.colab.kernel.invokeFunction("text_callback", [row, newValue], {});
}
function toggleCellContent(cell, row, col) {
  if (row === col) {
    var textCell = document.getElementById(`cell_${row}`);
    var diagMap = { 'H': ['E', 'yellow', 5], 'E': ['C', 'lime', 3], 'C': ['?', 'lightgray', 0], '?': ['H', 'red', 19],  };
    [cell.textContent, cell.style.backgroundColor, textCell.value] = diagMap[cell.textContent] || ['', ''];
    if (cell.textContent === "?" || cell.textContent === "H" || cell.textContent === "C"){
      var gridContainer = document.querySelector('.grid-container');
      var currentSize = (Math.sqrt(1 + 4 * gridContainer.childElementCount) - 1) / 2 - 1;
      for (var k = 0; k < currentSize; k++){
        if (k !== row){
          var a = document.getElementById(`cell_${k}_${row}`);
          var b = document.getElementById(`cell_${row}_${k}`);
          var diagCell = document.getElementById(`cell_${k}_${k}`);
          if (cell.textContent === "?" || diagCell.textContent === "?" ){
            a.textContent = b.textContent = '?';
            a.style.backgroundColor = b.style.backgroundColor = 'lightgray';
            a.style.opacity = b.style.opacity = 0.1;
          } else {
            a.textContent = b.textContent = '0';
            a.style.backgroundColor = b.style.backgroundColor = 'white';
            if (cell.textContent === "C" || diagCell.textContent === "C" ){
              a.style.opacity = b.style.opacity = 0.1;
            } else {
              a.style.opacity = b.style.opacity = 1.0;
            }
          }
        }
      }
    }
  } else {
    var diagRow = document.getElementById(`cell_${row}_${row}`);
    var diagCol = document.getElementById(`cell_${col}_${col}`);
    if (diagRow.textContent !== "?" && diagCol.textContent !== "?" && diagRow.textContent !== "C" && diagCol.textContent !== "C"){
      var symcell = document.getElementById(`cell_${col}_${row}`);
      var offdiagMap = { '0': ['1', 'lightblue'], '1': ['?', 'lightgray'], '?': ['0', 'white'] };
      [cell.textContent, cell.style.backgroundColor] = [symcell.textContent, symcell.style.backgroundColor] = offdiagMap[cell.textContent] || ['', ''];
    }
  }
  google.colab.kernel.invokeFunction("toggle_callback", [row, col], {});
}
function updateGridClickEvents() {
  var cells = document.querySelectorAll('.grid-item');
  cells.forEach(cell => {
    cell.addEventListener('click', () => {
      var [row, col] = cell.id.split('_').slice(1).map(Number);
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
  var positionInput = document.getElementById('position');
  var gridContainer = document.querySelector('.grid-container');
  var currentSize = (Math.sqrt(1 + 4 * gridContainer.childElementCount) - 1) / 2 - 1;
  var newSize = add ? currentSize + 1 : currentSize - 1;
  var position = parseInt(positionInput.value);
  if (position === -1) { position = currentSize; }

  var newGrid = '<div class="pos"></div>';
  for (var row = 0; row < newSize; row++) {
    newGrid += `<div class="pos">${row}</div>`;
  }
  newGrid += '<div class="pos"></div>';
  for (var row = 0; row < newSize; row++) {
    var oldRow = add ? (row < position ? row : row - 1) : (row < position ? row : row + 1);
    newGrid += `<div class="pos">${row}</div>`;
    for (var col = 0; col < newSize; col++) {
      var oldCol = add ? (col < position ? col : col - 1) : (col < position ? col : col + 1);
      var bgColor, content, opacity;
      if (add && (row === position || col === position)) {
        if (row === col) {
          bgColor = 'red';
          content = 'H';
          opacity = 1.0;
        } else {
          var cell = col === position ? document.getElementById(`cell_${oldRow}_${oldRow}`) : document.getElementById(`cell_${oldCol}_${oldCol}`);
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
        var cell = document.getElementById(`cell_${oldRow}_${oldCol}`);
        newGrid += createGridItem(row, col, cell.style.backgroundColor, cell.textContent, cell.style.opacity);
      }
    }
    if (add && row === position) {
      newGrid += createTextInput(row, 19);
    } else {
      var cell = document.getElementById(`cell_${oldRow}`);
      newGrid += createTextInput(row, cell.value);
    }
  }
  gridContainer.style.gridTemplateColumns = `repeat(${newSize + 2}, 30px)`;
  gridContainer.innerHTML = newGrid;
  updateGridClickEvents();
  positionInput.max = newSize;
  google.colab.kernel.invokeFunction("update_callback", [position, add], {});
}
updateGridClickEvents();