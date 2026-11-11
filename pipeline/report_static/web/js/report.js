$(function () {
    $(".bs-docs-sidenav>li>a").click(function () {
        maoTop(160)
    })

    $(".nav-stacked a").click(function () {
        maoTop(80)
    })

    $(".note-btn").click(function () {
        let note = $(this).next(".note")
        if (note.css("display") === "none") {
            $(note).show()
            let target_top = $(note).offset().top - 300;
            $("html,body").animate({scrollTop: target_top}, 800);
        } else {
            $(note).slideUp();
        }
    })


    //  $(".nav-stacked").hide();
    $(".bs-docs-sidenav > li").click(function () {
        // $(".bs-docs-sidenav > li").removeClass("active")
        $(".nav-stacked").hide();
        $(this).find(".nav-stacked").slideDown();
        $(".bs-docs-sidenav > li > i").removeClass("fa-caret-up").addClass(" fa-caret-down");
        $(this).find('i').removeClass("fa-caret-down").addClass("fa-caret-up");
    });

    $('body').scrollspy({target: '#myScrollspy'})

    $.each(tableJson, (i, v) => {
        initTable(v.id,v.info)
    })

    $.each(tableJsonNoSearch, (i, v) => {
        initTableNoSearch(v.id,v.info)
    })


})


function initTable( table,info) {

            let columns = []
            $.each(info.header, (i, v) => {
                columns.push({
                    field: `${i}`,
                    title: v.title,
                    sortable: true,
                    searchType: v.searchType
                })
            })

            let json = []

            let lines =  info.data.split("\n")
            $.each(lines, (i, v) => {
                let line = v.split("\t")
                let tableLine = {}
                $.each(line, (j, z) => {
                    tableLine[`${j}`] = z
                })
                json.push(tableLine)
            })

            $(`#${table}`).bootstrapTable({
                classes: "table table-hover table-border",
                pagination: "true",
                pageList: [10, 25, 50, 100, "All"],
                advancedSortable: "true",
                data: json,
                columns: columns

            })

}

function initTableNoSearch( table,info) {

    let columns = []
    $.each(info.header, (i, v) => {
        columns.push({
            field: `${i}`,
            title: v.title
        })
    })

    let json = []

    let lines =  info.data.split("\n")
    $.each(lines, (i, v) => {
        let line = v.split("\t")
        let tableLine = {}
        $.each(line, (j, z) => {
            tableLine[`${j}`] = z
        })
        json.push(tableLine)
    })

    $(`#${table}`).bootstrapTable({
        classes: "table table-hover table-border",
        data: json,
        columns: columns

    })

}

function maoTop(num) {
    window.onscroll = function () {
        let x = location.hash
        if (x) {
            $('html, body').animate({scrollTop: $(x).offset().top - num}, 0);
        }
    }
    setTimeout(function () {
        window.onscroll = null;
    }, 100)
}
