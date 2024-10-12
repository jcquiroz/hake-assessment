slide_buttons <- function(slide_id) {
  glue::glue('<p class="text-center buttons"><a class="btn btn-danger" target="_blank" href="{slide_id}.html"><i class="fa-solid fa-arrow-up-right-from-square"></i> Ver todos los diapositivas en una nueva ventana</a> <a class="btn btn-danger" target="_blank" href="{slide_id}.pdf" role="button"><i class="fa-solid fa-file-pdf"></i> Descargar PDF de todas las diapositivas</a></p>')
}
