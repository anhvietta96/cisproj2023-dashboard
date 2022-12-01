
var dynamicchart_bool = document.getElementById("dynamic-chart");
if(dynamicchart_bool){
  const data = JSON.parse(document.currentScript.nextElementSibling.textContent);
  console.log(data);
}