
var sidebarOpen = false
var sidebar = document.getElementById("sidebar")


function openSidebar() {
    if(!sidebarOpen)
        sidebar.classList.add("sidebar-responsive");
        sidebarOpen = true;
}

function closeSidebar() {
    if(sidebarOpen)
        sidebar.classList.remove("sidebar-responsive");
        sidebarOpen = false;
}

var histogramDom = document.getElementById('histogram');
var histogram = echarts.init(histogramDom,null);
var histogramOption = {
  legend: {
    data: ['Drugs'],
    right: '10%',
    top: '3%',
  },
  series: [
    {
      barGap: '0%',
      data: [
        [-1.775,0],
        [-1.325,0],
        [-0.875,3],
        [-0.425,2],
        [0.025,0],
        [0.475,0],
        [0.925,4],
        [1.375,3],
        [1.825,2],
        [2.275,1],
        [3.175,2],
        [3.625,3],
        [4.075,2],
        [4.525,1],
        [4.975,2],
        [5.425,0],
        [5.875,1],
        [6.325,1],
        [6.775,0],
      ],
      emphasis:{
        focus:'self',
        label: {
          show:true,
          color:'black',
          offset: [0,-100],
          fontSize: 15,
          backgroundColor: 'white',
          borderColor: 'black',
          borderWidth: 2,
          borderType: 'solid',
          padding:5,
          formatter: function(param){
            return `${param.data[0]},${param.data[1]}`;
          }
        }
      },
      type: 'bar',
    }
  ],
  xAxis: {
    type:'value',
    min: -2.5,
    max: 8.5,
    name: 'Attribute',
    nameLocation: 'start',
    nameTextStyle:{
      fontWeight:'bold',
      fontSize: 12,
      borderType: 'solid',},
    nameGap:30,},
  yAxis: {
    type:'value',
    min: 0,
    max: 5,
    name: 'Count',
    nameLocation: 'middle',
    nameTextStyle:{
      fontWeight:'bold',
      fontSize: 12,
      borderType: 'solid',},
    nameGap:50},
};

histogramOption && histogram.setOption(histogramOption);
console.log(histogramOption);

var chartDom = document.getElementById('scatter-chart');
var myChart = echarts.init(chartDom,null);
var option = {
  legend: {
    data: ['Drugs'],
    right: '10%',
    top: '3%',
  },
  series: [
    {
      data: [
        [1.6109, 474.587], 
        [1.1981, 285.343],
        [1.0482, 315.369],
        [1.3101, 180.159],
        [2.3753, 547.674],
        [5.7358, 314.469],
        [-0.7976, 324.341],
        [3.1869, 378.108],
        [-0.3513, 297.745],
        [4.1686, 277.411],
        [1.09718, 499.534],
        [4.4256, 375.871],
        [0.76592, 311.363],
        [3.809, 412.946],
        [0.9744, 312.391],
        [3.7357, 853.918],
        [1.5657, 193.246],
        [1.6344, 449.385],
        [4.8106, 319.88],
        [3.81298, 324.3993],
        [6.3136, 558.65],
        [-1.0397, 180.167],
        [3.9624, 311.469],
        [-0.4636, 366.084],
        [2.65458, 425.754],
        [3.1482, 388.895],
        [4.8944, 318.873],
      ],
      emphasis:{
        focus:'self',
        label: {
          show:true,
          color:'black',
          offset: [0,-100],
          fontSize: 15,
          backgroundColor: 'white',
          borderColor: 'black',
          borderWidth: 2,
          borderType: 'solid',
          padding:5,
          formatter: function(param){
            return `${param.data[0]},${param.data[1]}`;
          }
        }
      },
      symbolSize: 15,
      type: 'scatter',
    }
  ],
  xAxis: {
    type:'value',
    min: -3,
    max: 8,
    name: 'Attribute 1',
    nameLocation: 'start',
    nameTextStyle:{
      fontWeight:'bold',
      fontSize: 12,
      borderType: 'solid',},
    nameGap:30,},
  yAxis: {
    type:'value',
    min: -89,
    max: 1123,
    name: 'Attribute 2',
    nameLocation: 'middle',
    nameTextStyle:{
      fontWeight:'bold',
      fontSize: 12,
      borderType: 'solid',},
    nameGap:50},
};

option && myChart.setOption(option);

window.onresize = function() {
  myChart.resize();
  histogram.resize();
};