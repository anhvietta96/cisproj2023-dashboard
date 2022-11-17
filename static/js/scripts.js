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

var areaChartOptions = {
    series = [{
        name: 'Var 1',
        data: [31,40,28,51,42,109,200]
    }],
    chart: {
            height: 350,
            type: 'area',
            toolbar:{
                show: false,
            },
        },
        color: ["#4f35a1"],
        dataLabels: {
            enabled: false,
        },
        stroke:{
            curve:"smooth"
        },
        label:["X1","X2","X3","X4","X5","X6","X7"],
        markers:{
            size: 0
        },
        yaxis:[
            {
                title:{
                    text:'Value',
                },
            }
        ],
        tooltip: {
            shared: true,
            intersect: false,
        }
};

var = areaChart = new ApexCharts(document.querySelector("#area-chart"),areaChartOptions);
areaChart.render()