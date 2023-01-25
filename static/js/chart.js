function checkValidInt(input)
{
 return typeof(input)=="string" && !isNaN(input) && !isNaN(parseInt(input));
}

function checkValidFloat(input)
{
  return typeof(input)=="string" && !isNaN(input) && !isNaN(parseFloat(input));
}

function sanitizeFloatInput(input,positive_bool)
{
  if(typeof(input)!="string")
  {
    
    return 0;
  }
  if(isNaN(input) || isNaN(parseFloat(input)))
  {
    if(positive_bool)
    {
      return Infinity;
    }
    else
    {
      return -Infinity;
    }
  }
  return parseFloat(input);
}

function numerical_filter(data,pos,lower_bound,upper_bound)
{
  var filtered_data = [];
  for(const i of Object.keys(data))
  {
    var sub_arr = [];
    for(let j = 0;j<data[i].length;j++)
    {
      if(data[i][j][pos]>lower_bound && data[i][j][pos]<upper_bound)
      {
        sub_arr.push(data[i][j]);
      }
    }
    filtered_data.push(sub_arr);
  }
  return filtered_data;
}

function lexicographic_filter(data,pos,substring)
{
  var filtered_data = [];
  for(const i of Object.keys(data))
  {
    var sub_arr = [];
    for(let j = 0;j<data[i].length;j++)
    {
      if(data[i][j][pos].includes(substring))
      {
        sub_arr.push(data[i][j]);
      }
    }
    filtered_data.push(sub_arr);
  }
  return filtered_data;
}

function create_new_input_field(type,id_arr)
{
  var sub_paragraph = document.createElement('a');
  var a_id;
  var input_1_id;
  var input_2_id;
  if(id_arr.length == 1)
  {
    a_id = `dynamic_input_${id_arr[0]}_`;
    input_1_id = `filter_input_${id_arr[0]}_0`;
    input_2_id = `filter_input_${id_arr[0]}_1`;
  }
  if(id_arr.length == 2)
  {
    a_id = `subparagraph-${id_arr[0]}-${id_arr[1]}-`;
    input_1_id = `subfilter-input-${id_arr[0]}-${id_arr[1]}-0`;
    input_2_id = `subfilter-input-${id_arr[0]}-${id_arr[1]}-1`;
  }
  if(type == 0)
  {
    sub_paragraph.id = a_id+type;
    sub_paragraph.appendChild(document.createTextNode(' contains '));
    var input_1 = document.createElement('input');
    input_1.id = input_1_id;
    sub_paragraph.appendChild(input_1);
  }
  if(type == 1)
  {
    sub_paragraph.id = a_id+type;
    sub_paragraph.appendChild(document.createTextNode(' ranges from '));
    var input_1 = document.createElement('input');
    input_1.id = input_1_id;
    sub_paragraph.appendChild(input_1);
    sub_paragraph.appendChild(document.createTextNode(' to '));
    var input_2 = document.createElement('input');
    input_2.id = input_2_id;
    sub_paragraph.appendChild(input_2);
  }
  return sub_paragraph;
}

function filter_multiple_attr(collected_filter,curr_data,threshold)
{
  var filtered = {}
  for(const group of Object.keys(curr_data))
  {
    filtered[group] = [];
    for(const i of Object.keys(curr_data[group]))
    {
      var data_point = curr_data[group][i];
      var count = 0;
      for(let j = 0; j < collected_filter.length; j++)
      {
        var attr = parseInt(collected_filter[j][0]);
        if(collected_filter[j].length == 2)
        {
          if(data_point[attr].includes(collected_filter[j][1]))
          {
            count++;
          }
        }
        else
        {
          if(collected_filter[j][1] < data_point[attr] && data_point[attr] < collected_filter[j][2])
          {
            count++;
          }
        }
      }
      if(count >= threshold)
      {
        filtered[group].push(data_point);
      }
    }
  }
  return filtered;
}

function get_selected_inchikey(option,chart_type)
{
  if([1,2].includes(chart_type))
  {
    var position = 2;
  }
  else if([3].includes(chart_type))
  {
    var position = 1;
  }
  var inchikey_array = [];
  for(const group of Object.keys(option['series']))
  {
    for(let i = 0; i < option['series'][group]['data'].length; i++)
    {
      inchikey_array.push(option['series'][group]['data'][i][position]);
    }
  }
  return inchikey_array;
}

function get_new_max_vals(data)
{
  var xy_min_max = [null,null,null,null];
  let x_values = [];
  let y_values = [];
  let x_col = [];
  let y_col = [];
  for(const i of Object.keys(data))
  {
    x_col = data[i].map(function(value,index){return value[0];});
    x_values.push.apply(x_values,x_col);
    y_col = data[i].map(function(value,index){return value[1];});
    y_values.push.apply(y_values,y_col);
  }
  xy_min_max[0] = Math.min.apply(null,x_values);
  xy_min_max[1] = Math.max.apply(null,x_values);
  xy_min_max[2] = Math.min.apply(null,y_values);
  xy_min_max[3] = Math.max.apply(null,y_values);
  x_dist = (xy_min_max[1] - xy_min_max[0])/5;
  y_dist = (xy_min_max[3] - xy_min_max[2])/2.5;

  xy_min_max[0] = Math.round(xy_min_max[0]-x_dist)-1;
  xy_min_max[1] = Math.round(xy_min_max[1]+x_dist)+1;
  xy_min_max[2] = Math.round(xy_min_max[2]-y_dist)-1;
  xy_min_max[3] = Math.round(xy_min_max[3]+y_dist)+1;

  return xy_min_max;
}

function extract_column(data,position)
{
  var separated_col = [];
  var concat_col = [];
  for(const i of Object.keys(data))
  {
    separated_col.push(null)
    separated_col[i]=data[i].map(function(value,index){return value[position];});
    concat_col.push.apply(concat_col,separated_col[i]);
  }
  return {separated_col:separated_col,concat_col:concat_col};
}

function get_scatter_option(data,image,legend,tooltip_col)
{
  var get_rich_text = function(chart_data,image)
  {
    var rich_text = [];
    var curr_rich_text = {Attr:{width:40,align:'left'},Val:{width:60,align:'right'}};
    for(const i of Object.keys(chart_data))
    {
      if(data_size < data_size_threshold)
      {
        for(const j of Object.keys(chart_data[i]))
        {
          var inchikey = chart_data[i][j][2];
          var img_url = image[inchikey];
          var inchikey_reformatted = inchikey.replace(/-/g,'');
          curr_rich_text[`${inchikey_reformatted}`]={height:150,backgroundColor:{image:img_url}}
        }
      }
      rich_text.push(curr_rich_text);
      curr_rich_text = {Attr:{width:40,align:'left'},Val:{width:60,align:'right'}};
    }
    return rich_text;
  }
  var get_series = function(chart_data,tooltip_col,rich_text)
  {
    var series_val = [];
    var choice_Min;
    var choice_diff;
    if(chart_type==2)
    {
      var choice_col = extract_column(chart_data,choice_param).concat_col;
      var choice_Max = Math.round(Math.max.apply(null,choice_col));
      choice_Min = Math.round(Math.min.apply(null,choice_col));
      choice_diff = choice_Max-choice_Min;   
    }
    for(const i of Object.keys(chart_data))
    {
      var series = {
        name: ChartOptions['legend'][i],
        data: chart_data[i],
        type: 'scatter',
        symbolSize: function(param){
          if(chart_type==1){return 15;}
          else{
          var minSize = 10;
          var maxSize = 40;
          var Sizediff = maxSize-minSize;
          var symbolSz = (param[choice_param]-choice_Min)*Sizediff/choice_diff+minSize;
          return symbolSz;}
        },
        emphasis: {
          focus: 'self',
          label: {
            show: true,
            formatter: function (param) {
              var return_list;
              var inchikey = param.data[2].replace(/-/g,'');
              if(data_size<data_size_threshold)
                {return_list = [`{${inchikey}|}`];}
              else
                {return_list = [];}
              for(const i of Object.keys(tooltip_col))
              {
                if(!lexicographic_position.includes(parseInt(i)))
                {
                  return_list.push(`{Attr|${tooltip_col[i]}}{Val|${param.data[i]}}`);
                }
              }
              var rtext = return_list.join('\n');
              return rtext;
            },
            position: 'top',
            color: '#000',
            backgroundColor: '#fff',
            fontWeight: 'bold',
            padding: [3,4],
            width: 300,
            height: 100,
            fontSize: 10,
            align: 'center',
            verticalAlign: 'bottom',
            rich:rich_text[i],
          }
        },
      }
      series_val.push(series);
    }
    return series_val;
  }
  
  var rich_text = get_rich_text(data,image);
  var series_val = get_series(data,tooltip_col,rich_text);
  var xy_min_max = get_new_max_vals(data);

  var option={
      legend: {
        right: '10%',
        top: '3%',
        data: legend,
      },
      grid:{
        show: true,
      },
      xAxis: {
        type:'value',
        min: xy_min_max[0],
        max: xy_min_max[1],
        name: tooltip_col[0],
        nameLocation: 'start',
        nameTextStyle:{
          fontWeight:'bold',
          fontSize: 12,
          borderType: 'solid',},
        nameGap:30,
        axisTick:{
          show: true,
          interval: 'auto',
        },},
      yAxis: {
        type:'value',
        min: xy_min_max[2],
        max: xy_min_max[3],
        name: tooltip_col[1],
        nameLocation: 'middle',
        nameTextStyle:{
          fontWeight:'bold',
          fontSize: 12,
          borderType: 'solid',},
        nameGap:50},
      series: series_val,
    };

  return option;
}

function process_histogram_data(data, num_of_bar)
{
  var extracted_col = extract_column(data,0);
  x_col = extracted_col.separated_col;
  x_values = extracted_col.concat_col;
  var x_min = Math.floor(Math.min.apply(null,x_values));
  var x_max = Math.floor(Math.max.apply(null,x_values))+1;
    
  let interval = (x_max-x_min)/num_of_bar;
  let threshold_arr = [];
  for(let i = 0; i <= num_of_bar; i++)
  {
    threshold_arr.push(x_min + i*interval);
  }
  let x_axis_val = [];
  for(let i = 0; i < num_of_bar; i++)
  {
    x_axis_val.push((threshold_arr[i]+threshold_arr[i+1])/2);
  }


  let count_arr = [];
  for(const i of Object.keys(data))
  {
    count_arr.push([]);
    for(let j = 0; j < num_of_bar; j++)
    {
      count_arr[i].push(0);
    }
  }
  for(const k of Object.keys(data))
  {
    for(let i = 0; i < x_col[k].length; i++)
    {
      for(let j = 0; j < num_of_bar; j++)
      {
        if(x_col[k][i] >= threshold_arr[j] && x_col[k][i] < threshold_arr[j+1])
       {
          count_arr[k][j]++;
        }
      }
    }
  }

  var aggregated_count = [];
  for(const i of Object.keys(data))
  {
    aggregated_count.push.apply(aggregated_count,count_arr[i]);
  }
  var y_max = Math.max.apply(null,aggregated_count);


  let processed_data = [];
  for(const i of Object.keys(data))
  {
    processed_data.push([]);
    for(let j = 0; j < num_of_bar; j++)
    {
      processed_data[i].push([x_axis_val[j],count_arr[i][j],Number(threshold_arr[j].toFixed(2)),Number((threshold_arr[j+1].toFixed(2)))]);
    }
  }

  return {processed_data:processed_data,x_min:x_min,x_max:x_max,y_max:y_max};
}



function get_histogram_option(data,legend,tooltip_col,num_of_bar)
{
  var histogram_data = process_histogram_data(data,num_of_bar);
  var series_val= []
  for(const i of Object.keys(data))
  {
    var series = {
      name: ChartOptions['legend'][i],
      data: histogram_data.processed_data[i],
      type: 'bar',
      barGap: '0%',
      emphasis:{
        focus: 'self',
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
            return [`Range [${param.data[2]},${param.data[3]}]`,`Count ${param.data[1]}`].join('\n');
          }
        }
      }
    }
    series_val.push(series);
  }

  option={
    legend: {
      right: '10%',
      top: '3%',
      data: ChartOptions['legend'],
    },
    xAxis: {
      type:'value',
      min: histogram_data.x_min,
      max: histogram_data.x_max,
      name: ChartOptions['header'][0],
      nameLocation: 'start',
      nameTextStyle:{
        fontWeight:'bold',
        fontSize: 12,
        borderType: 'solid',},
      nameGap:30,},
    yAxis: {
      type:'value',
      min: 0,
      max: histogram_data.y_max,
      name: 'Count',
      nameLocation: 'middle',
      nameTextStyle:{
        fontWeight:'bold',
        fontSize: 12,
        borderType: 'solid',},
      nameGap:50},
    series: series_val,
  };
  return option;
}

var ChartOptions;
var Chart;
var option;
var chart_type;
var lexicographic_position;
var curr_num_of_bar = 20;
var legend;
var tooltip_col;
var curr_data;
var choice_param = 5;
var data_size;
var data_size_threshold = 1000;

var chartDom = document.getElementById('dynamic-chart');
if(chartDom) {
  Chart = echarts.init(chartDom,null);
  ChartOptions = JSON.parse(document.currentScript.nextElementSibling.textContent);
  tooltip_col = ChartOptions['header'];
  curr_data = $.extend(true,{},ChartOptions['data']);
  var image = ChartOptions['image'];
  legend = ChartOptions['legend'];
  chart_type = ChartOptions['type'];
  lexicographic_position = ChartOptions['lexicographic_position'];
  data_size = ChartOptions['size'];
  
  if([1,2].includes(chart_type))
  {
    option = get_scatter_option(curr_data,image,legend,tooltip_col);
  }
  else if([3].includes(chart_type))
  {
    option = get_histogram_option(curr_data,legend,tooltip_col,curr_num_of_bar);
  }
  option && Chart.setOption(option);

  var csv_input = document.getElementById('export-csv-val');
  csv_input.value = JSON.stringify(get_selected_inchikey(option,chart_type));

  var sdf_input = document.getElementById('export-sdf-val');
  sdf_input.value = csv_input.value;
}

$(document).ready(function(){
  var num_filter = 0;
  var available_filter = [];
  var available_custom_filter = [];
  var custom_filter_tracker = {};
  $('#add-filter').click(function(){
    var index;
    for(let i = 0; i <= num_filter;i++)
    {
      if(!available_filter.includes(i))
      {
        index = i;
        break;
      }
    }
    var new_div = document.createElement('div');
    new_div.id='filter_div_'+index;
    var new_select = document.createElement('select');
    new_select.id='filter_'+index;
    var header = ChartOptions['header'];
    var html = '';
    for(let i = 0;i<header.length;i++)
    {
      html += `<option value=${i}>${header[i]}</option>`;
    }
    var html_a=`<a id='dynamic_input_${index}_1'>`;
    html_a += ` ranges from <input type=text id=filter_input_${index}_0></input> to <input type=text id=filter_input_${index}_1></input>`;

    var button = document.createElement('button');
    button.id='remove_'+index;
    button.classList.add('filter-btn');
    var close_span = document.createElement('span');
    close_span.classList.add('material-icons-outlined');
    close_span.appendChild(document.createTextNode('close'));
    button.appendChild(close_span);
    const element = document.getElementById('filters');
    const elem_button = document.getElementById('add-filter');
    new_div.appendChild(button);
    new_div.appendChild(new_select);
    element.insertBefore(new_div,elem_button);
    $('#filter_'+index).append(html);
    $('#filter_div_'+index).append(html_a);
    num_filter++;
    available_filter.push(index);
    return;
  });
  $('#add-custom-filter').click(function(){
    var index;
    for(let i = 0; i <= available_custom_filter.length ;i++)
    {
      if(!available_custom_filter.includes(i))
      {
        index = i;
        break;
      }
    }
    available_custom_filter.push(index);
    custom_filter_tracker[index] = [];

    var new_div = document.createElement('div');
    new_div.id = `custom_filter_div_${index}`;
    var button = document.createElement('button');
    button.id = `remove_custom_${index}`;
    button.classList.add('sub-filter-btn');
    var close_span = document.createElement('span');
    close_span.classList.add('material-icons-outlined');
    close_span.appendChild(document.createTextNode('close'));
    button.appendChild(close_span);
    var main_div_text = document.createTextNode(' Data points must satisfy ');
    var custom_filter_input = document.createElement('input');
    custom_filter_input.id = `custom_filter_input_${index}`;
    custom_filter_input.type = 'number';
    custom_filter_input.value = 0;
    custom_filter_input.size = 2;
    custom_filter_input.min = 0;
    custom_filter_input.max = 0;
    var total_filter_text = document.createElement('a');
    total_filter_text.id = `total_filter_${index}`;
    total_filter_text.appendChild(document.createTextNode(' out of 0 filters\n'));
    var add_filter_btn = document.createElement('button');
    add_filter_btn.id = `add_filter_rule_${index}`;
    add_filter_btn.classList.add('custom_set_btn');
    add_filter_btn.appendChild(document.createTextNode('Add Filter in Custom Set'));
    
    new_div.appendChild(button);
    new_div.appendChild(main_div_text);
    new_div.appendChild(custom_filter_input);
    new_div.appendChild(total_filter_text);
    new_div.appendChild(add_filter_btn);

    const element = document.getElementById('filters');
    const elem_button = document.getElementById('add-filter');
    element.insertBefore(new_div,elem_button);
  });
  $('#filters').on('click','button.custom_set_btn', function(){
    var id = this.id;
    var index = id.split("_")[3];

    var custom_filter_input = document.getElementById(`custom_filter_input_${index}`);
    custom_filter_input.max++;
    var curr_max = custom_filter_input.max;
    document.getElementById(`total_filter_${index}`).remove();
    var total_filter_text = document.createElement('a');
    total_filter_text.id = `total_filter_${index}`;
    total_filter_text.appendChild(document.createTextNode(` out of ${curr_max} filters\n`));
    document.getElementById(`custom_filter_div_${index}`).insertBefore(total_filter_text,document.getElementById(`add_filter_rule_${index}`));

    var new_rule_index;
    for(let i = 0; i <= custom_filter_tracker[index].length; i++)
    {
      if(!custom_filter_tracker[index].includes(i))
      {
        new_rule_index = i;
      }
    }
    custom_filter_tracker[index].push(new_rule_index);


    var sub_div = document.createElement('div');
    sub_div.style.paddingLeft = '30px';
    sub_div.id = `rule_div_${index}_${new_rule_index}`;
    var sub_delete_button = document.createElement('button');
    sub_delete_button.id = `rule_${index}_${new_rule_index}`;
    sub_delete_button.classList.add('sub-delete-btn');
    var close_span = document.createElement('span');
    close_span.classList.add('material-icons-outlined');
    close_span.appendChild(document.createTextNode('close'));
    sub_delete_button.appendChild(close_span);
    var sub_select = document.createElement('select');
    sub_select.id = `sub_select_${index}_${new_rule_index}`;
    for(const i of Object.keys(tooltip_col))
    {
      var sub_option = document.createElement('option');
      sub_option.value = i;
      var text_node = document.createTextNode(tooltip_col[i]);
      sub_option.appendChild(text_node);
      sub_select.appendChild(sub_option);
    }
    var sub_paragraph = document.createElement('a');
    sub_paragraph.id = `subparagraph-${index}-${new_rule_index}-1`;
    sub_paragraph.appendChild(document.createTextNode(' ranges from '));
    var input_1 = document.createElement('input');
    input_1.id = `subfilter-input-${index}-${new_rule_index}-0`;
    sub_paragraph.appendChild(input_1);
    sub_paragraph.appendChild(document.createTextNode(' to '));
    var input_2 = document.createElement('input');
    input_2.id = `subfilter-input-${index}-${new_rule_index}-1`;
    sub_paragraph.appendChild(input_2);

    sub_div.appendChild(document.createTextNode('\t'));
    sub_div.appendChild(sub_delete_button);
    sub_div.appendChild(sub_select);
    sub_div.appendChild(sub_paragraph);

    var main_div = document.getElementById(`custom_filter_div_${index}`);
    main_div.appendChild(sub_div);
  });
  $('#filters').on('change','select',function(events){
    var id = this.id;
    var split_id = id.split("_");
    if(split_id[0] == "filter")
    {
      var index = split_id[1];
      if(document.getElementById(`dynamic_input_${index}_1`) && lexicographic_position.includes(parseInt(this.value)))
      {
        document.getElementById(`dynamic_input_${index}_1`).remove();
        var new_field = create_new_input_field(0,[index]);
        document.getElementById(`filter_div_${index}`).appendChild(new_field);
      }
      if(document.getElementById(`dynamic_input_${index}_0`) && !lexicographic_position.includes(parseInt(this.value)))
      {
        document.getElementById(`dynamic_input_${index}_0`).remove();
        var new_field = create_new_input_field(1,[index]);
        document.getElementById(`filter_div_${index}`).appendChild(new_field);
      }
    }
    else
    {
      var set_index = split_id[2];
      var rule_index = split_id[3];
      if(document.getElementById(`subparagraph-${set_index}-${rule_index}-1`) && lexicographic_position.includes(parseInt(this.value)))
      {
        document.getElementById(`subparagraph-${set_index}-${rule_index}-1`).remove();
        var new_field = create_new_input_field(0,[set_index,rule_index]);
        document.getElementById(`rule_div_${set_index}_${rule_index}`).appendChild(new_field);
      }
      if(document.getElementById(`subparagraph-${set_index}-${rule_index}-0`) && !lexicographic_position.includes(parseInt(this.value)))
      {
        document.getElementById(`subparagraph-${set_index}-${rule_index}-0`).remove();
        var new_field = create_new_input_field(1,[set_index,rule_index]);
        document.getElementById(`rule_div_${set_index}_${rule_index}`).appendChild(new_field);
      }
    }
  });
  $('#filters').on('click','button.filter-btn',function(events){
    var id = this.id;
    var split_id = id.split("_");
    var deleteindex = split_id[1];
    $('#filter_div_'+deleteindex).remove();
    num_filter--;
    var new_arr = [];
    for(let i = 0;i<available_filter.length;i++)
    {
      if(available_filter[i]!=deleteindex)
      {
        new_arr.push(available_filter[i]);
      }
    }
    available_filter=new_arr;
    return;
  });
  $('#filters').on('click','button.sub-filter-btn',function(){
    var id = this.id;
    var split_id = id.split("_");
    var deleteindex = split_id[2];
    $('#custom_filter_div_'+deleteindex).remove();

    var arr_index = available_custom_filter.indexOf(parseInt(deleteindex));
    if(arr_index !== -1)
    {
      available_custom_filter.splice(arr_index,1);
      delete custom_filter_tracker[arr_index];
    }
    return;
  });
  $('#filters').on('click','button.sub-delete-btn',function(){
    var id = this.id;
    var split_id = id.split("_");
    var set_index = split_id[1];
    var rule_index = split_id[2];
    $(`#rule_div_${set_index}_${rule_index}`).remove();

    var arr_index = custom_filter_tracker[set_index].indexOf(parseInt(rule_index));

    if(arr_index !== -1)
    {
      custom_filter_tracker[set_index].splice(arr_index,1);
    }
    return;
  });
  $('#apply-filter').click(function(){
    curr_data = $.extend(true,{},ChartOptions['data']);
    if(available_filter.length!=0)
    {
      var filter = null;
      var filter_val = [null,null];
      for(const slot of available_filter)
      {
        filter = document.getElementById('filter_'+slot);
        if(!lexicographic_position.includes[filter.value])
        {
          filter_val[0] = document.getElementById(`filter_input_${slot}_0`).value;
          filter_val[1] = document.getElementById(`filter_input_${slot}_1`).value;
          curr_data = numerical_filter(curr_data,filter.value,sanitizeFloatInput(filter_val[0],false),sanitizeFloatInput(filter_val[1],true));
        }
        else
        {
          filter_val[0] = document.getElementById(`filter_input_${slot}_0`).value;
          curr_data = lexicographic_filter(curr_data,filter.value,filter_val[0]);
        }
      }
    }
    if(available_custom_filter.length != 0)
    {
      for(const set_index of Object.keys(custom_filter_tracker))
      {
        var threshold = document.getElementById(`custom_filter_input_${set_index}`).value;
        var collected_filter = [];
        for(let i = 0; i < custom_filter_tracker[set_index].length; i++)
        {
          var rule_index = custom_filter_tracker[set_index][i];
          var attribute = document.getElementById(`sub_select_${set_index}_${rule_index}`).value;
          if(lexicographic_position.includes(parseInt(attribute)))
          {
            let val_1 = document.getElementById(`subfilter-input-${set_index}-${rule_index}-0`).value;
            collected_filter.push([attribute,val_1]);
          }
          else
          {
            let val_1 = document.getElementById(`subfilter-input-${set_index}-${rule_index}-0`).value;
            let val_2 = document.getElementById(`subfilter-input-${set_index}-${rule_index}-1`).value;
            collected_filter.push([attribute,sanitizeFloatInput(val_1,false),sanitizeFloatInput(val_2,true)]);
          }
        }
        curr_data = filter_multiple_attr(collected_filter,curr_data,threshold);
      }
    }
    if(chart_type == 3)
    {
      var histogram_data = process_histogram_data(curr_data,curr_num_of_bar);
      for(const i of Object.keys(curr_data))
      {
        option['series'][i]['data'] = histogram_data.processed_data[i];
      }
      option['xAxis']['min']=histogram_data.x_min;
      option['xAxis']['max']=histogram_data.x_max;
      option['yAxis']['max']=histogram_data.y_max;
    }
    else if([1,2].includes(chart_type))
    {
      for(const i of Object.keys(curr_data))
      {
        option['series'][i]['data'] = curr_data[i];
      }
      var xy_min_max = get_new_max_vals(curr_data);
      option['xAxis']['min'] = xy_min_max[0];
      option['xAxis']['max'] = xy_min_max[1];
      option['yAxis']['min'] = xy_min_max[2];
      option['yAxis']['max'] = xy_min_max[3];
    }
    
    option && Chart.setOption(option);
    
    var csv_input = document.getElementById('export-csv-val');
    csv_input.value = JSON.stringify(get_selected_inchikey(option));

    var sdf_input = document.getElementById('export-sdf-val');
    sdf_input.value = csv_input.value;
  });
  $('#apply-customization').click(function(){
    var all_options = document.getElementById('customize-options');
    var all_custom_select = all_options.querySelectorAll('select');
    for(var i = 0; i < all_custom_select.length; i++)
    {
      var select_field = all_custom_select[i];
      if(select_field.id == 'new-param-choice')
      {
        if(select_field.value != -1)
        {
          chart_type = 2;
          choice_param = select_field.value-1;
          option = get_scatter_option(curr_data,image,legend,tooltip_col);
        }
        else
        {
          chart_type = 1;
          option = get_scatter_option(curr_data,image,legend,tooltip_col);
        }
      }
    }
    var all_custom_input = all_options.querySelectorAll('input');
    for(var i = 0; i < all_custom_input.length; i++)
    {
      var input_field = all_custom_input[i];
      if(input_field.id == 'new-min-x' && checkValidFloat(input_field.value))
      {
        option['xAxis']['min'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'new-max-x' && checkValidFloat(input_field.value))
      {
        option['xAxis']['max'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'new-min-y' && checkValidFloat(input_field.value))
      {
        option['yAxis']['min'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'new-max-y' && checkValidFloat(input_field.value))
      {
        option['yAxis']['max'] = parseFloat(input_field.value);
      }
      if(input_field.id == 'bin-num'  && checkValidInt(input_field.value))
      {
        curr_num_of_bar = parseInt(input_field.value);
        var histogram_data = process_histogram_data(curr_data,curr_num_of_bar);
        for(const i of Object.keys(curr_data))
        {
          option['series'][i]['data'] = histogram_data.processed_data[i];
        }
        option['xAxis']['min']=histogram_data.x_min;
        option['xAxis']['max']=histogram_data.x_max;
        option['yAxis']['max']=histogram_data.y_max;
      }
    }
    option && Chart.setOption(option);
  });
  $('#export-png').click(function(){
    var canvas_collection = document.getElementsByTagName('canvas');
    var dataURL = canvas_collection[0].toDataURL("image/png");
    var newTab = window.open('about:blank','image from canvas');
    newTab.document.write("<img src='"+dataURL+"' alt='from canvas'/>");
  });
});

window.onresize = function() {
  Chart.resize();
};